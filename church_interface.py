%%writefile {save_dir}church_interface.py


# -*- coding: utf-8 -*-
###
# MolFS Church Algorithm device
# Implementation of Church's DNA Storage Algorithm for MolFS

import os
import sys
import tempfile
import csv
import shutil
from church_algorithm_encoder import encode_file_to_dna
from church_algorithm_decoder import decode_oligos_from_csv, decode_oligos_from_fastq, BASE_TO_BIT

class MolFSDev:
    """
    Church Algorithm interface for MolFS
    
    This interface implements the Church DNA storage algorithm
    for use with the Molecular File System (MolFS).
    """
    
    # Default universal primers (can be overridden)
    DEFAULT_FORWARD_PRIMER = "CTACACGACGCTCTTCCGATCT"
    DEFAULT_REVERSE_PRIMER = "AGATCGGAAGAGCGGTTCAGCA"
    
    # Whether to use FastQ format (for future implementation)
    UseFastQ = False
    
    def __init__(self):
        """
        Initialize the Church Algorithm interface for MolFS
        """
        self.DoubleStep = False
        self.ValidDecode = False
        
        self.encodeParam = 0
        self.decodeParam = 0
        
        self.Address = 0
        self.Block = 0
        self.Pool = 0
        
        self.Redundancy = 0
        self.MaxIterations = 0
        
        # Block size configuration (in bytes)
        self.block_size_bytes = 5 * 1024  # Default: 5KB blocks
        
        # Create temporary directory for intermediate files
        self.temp_dir = tempfile.mkdtemp(prefix="church_molfs_")
        
        # Initialize primer registry
        # This allows external code to register primers for specific Pool/Block combinations
        self.primer_registry = {}
    
    def __del__(self):
        """
        Clean up temporary directory when the object is destroyed
        """
        try:
            if os.path.exists(self.temp_dir):
                shutil.rmtree(self.temp_dir)
        except:
            pass
    
    def set_block_size(self, size_bytes):
        """
        Set the block size for encoding
        
        Args:
            size_bytes (int): Block size in bytes
        """
        self.block_size_bytes = size_bytes
    
    def register_primers(self, pool, block, forward_primer, reverse_primer):
        """
        Register primers for a specific Pool/Block combination
        
        Args:
            pool (int): Pool number
            block (int): Block number
            forward_primer (str): Forward primer sequence
            reverse_primer (str): Reverse primer sequence
        """
        key = (pool, block)
        self.primer_registry[key] = (forward_primer, reverse_primer)
    
    def get_primers(self, pool, block):
        """
        Get primers for a specific Pool/Block combination
        
        Args:
            pool (int): Pool number
            block (int): Block number
            
        Returns:
            tuple: (forward_primer, reverse_primer)
        """
        key = (pool, block)
        
        # Return registered primers if available
        if key in self.primer_registry:
            return self.primer_registry[key]
        
        # Otherwise, return default primers
        return self.DEFAULT_FORWARD_PRIMER, self.DEFAULT_REVERSE_PRIMER
    
    def encode_block(self, in_file, out_file, block_index, start_offset=0, block_length=None):
        """
        Encode a specific block from a file using the current Pool/Block settings
        
        This allows encoding selective blocks from a file and placing them in
        different pools with different primers.
        
        Args:
            in_file (str): Path to the input binary file
            out_file (str): Path to save the DNA sequences
            block_index (int): Index of this block (for reconstruction)
            start_offset (int, optional): Byte offset where this block starts in the file
            block_length (int, optional): Length of this block in bytes (defaults to block_size_bytes)
            
        Returns:
            tuple: (Success status, Block information)
        """
        # Use the configured block size if not specified
        if block_length is None:
            block_length = self.block_size_bytes
            
        # Create a temporary file for this block
        temp_block_file = os.path.join(self.temp_dir, f"temp_block_{block_index}.bin")
        
        try:
            # Extract this block from the file
            with open(in_file, 'rb') as f:
                f.seek(start_offset)
                block_data = f.read(block_length)
                
            # Write to temporary file
            with open(temp_block_file, 'wb') as f:
                f.write(block_data)
                
            # Set the current block index
            self.Block = block_index
            
            # Get primers for the current Pool/Block
            forward_primer, reverse_primer = self.get_primers(self.Pool, self.Block)
            
            # Encode this block
            num_oligos, file_size, num_blocks = encode_file_to_dna(
                temp_block_file,
                out_file,
                forward_primer=forward_primer,
                reverse_primer=reverse_primer,
                block_size_bytes=block_length  # Use the actual block length
            )
            
            print(f"Encoded block {block_index} ({len(block_data)} bytes) into {num_oligos} DNA oligos")
            print(f"Block assigned to Pool {self.Pool}")
            
            return True, {
                "block_index": block_index,
                "pool": self.Pool,
                "block_size": len(block_data),
                "num_oligos": num_oligos
            }
            
        except Exception as e:
            print(f"Error encoding block {block_index}: {e}")
            return False, {}
        finally:
            # Clean up temporary file
            if os.path.exists(temp_block_file):
                os.remove(temp_block_file)
    
    def encode(self, in_file, out_file):
        """
        Encode a binary file into DNA sequences (legacy method)
        
        This maintains compatibility with the original interface but delegates to encode_block
        
        Args:
            in_file (str): Path to the input binary file
            out_file (str): Path to save the DNA sequences
            
        Returns:
            tuple: (Success status, Number of blocks)
        """
        # Simply delegates to encode_block with block index 0
        success, block_info = self.encode_block(in_file, out_file, 0)
        
        if success:
            return True, 1
        else:
            return False, 0
    
    def encode_file(self, in_file, out_dir, distribution_strategy=None):
        """
        Encode an entire file across multiple pools using a distribution strategy
        
        Args:
            in_file (str): Path to the input binary file
            out_dir (str): Directory to save the DNA sequence files
            distribution_strategy (callable, optional): Function that determines which pools a block goes to
                Function signature: distribution_strategy(block_index, total_blocks) -> list of pool indices
                
        Returns:
            dict: Mapping of (pool, block) to output file and block info
        """
        # Create output directory if it doesn't exist
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
            
        # Get file size
        file_size = os.path.getsize(in_file)
        
        # Calculate number of blocks
        total_blocks = (file_size + self.block_size_bytes - 1) // self.block_size_bytes
        
        # Default distribution strategy: 1 block per pool
        if distribution_strategy is None:
            distribution_strategy = lambda block_idx, total: [block_idx % 3 + 1]  # Default to pools 1-3
            
        # Track where each block is stored
        block_distribution = {}
        
        # Process each block
        for block_idx in range(total_blocks):
            # Calculate block offset and length
            start_offset = block_idx * self.block_size_bytes
            block_length = min(self.block_size_bytes, file_size - start_offset)
            
            # Determine which pools this block should go to
            target_pools = distribution_strategy(block_idx, total_blocks)
            
            for pool_idx in target_pools:
                # Set the current pool
                self.Pool = pool_idx
                
                # Create output filename
                out_file = os.path.join(out_dir, f"pool{pool_idx}_block{block_idx}.csv")
                
                # Encode this block
                success, block_info = self.encode_block(
                    in_file, 
                    out_file, 
                    block_idx,
                    start_offset,
                    block_length
                )
                
                if success:
                    block_distribution[(pool_idx, block_idx)] = {
                        "file": out_file,
                        "info": block_info
                    }
                    
        print(f"Encoded {file_size} bytes into {total_blocks} blocks across {len(set(p for p, _ in block_distribution.keys()))} pools")
        return block_distribution
    
    def decode(self, in_file, out_file):
        """
        Decode DNA sequences back to the original binary file
        
        Args:
            in_file (str): Path to the input file with DNA sequences (FASTQ or CSV)
            out_file (str): Path to save the decoded binary file
            
        Returns:
            tuple: (Success status, Block information)
        """
        self.ValidDecode = False
        
        # Get primers for the current Pool/Block
        forward_primer, reverse_primer = self.get_primers(self.Pool, self.Block)
        
        # Check file extension to determine how to handle it
        file_extension = os.path.splitext(in_file)[1].lower()
        
        try:
            if file_extension in ['.fastq', '.fq']:
                # Handle FASTQ file (raw sequencing data without consensus/error correction)
                decoded_blocks, block_info = decode_oligos_from_fastq(
                    in_file,
                    forward_primer=forward_primer,
                    reverse_primer=reverse_primer,
                    output_file=out_file
                )
            else:
                # Handle CSV file (processed data)
                decoded_blocks, block_info = decode_oligos_from_csv(
                    in_file,
                    forward_primer=forward_primer,
                    reverse_primer=reverse_primer,
                    output_file=out_file
                )
            
            if decoded_blocks:
                self.ValidDecode = True
                print(f"Decoded {len(decoded_blocks)} data blocks")
                print(f"Number of blocks: {block_info['block_count']}")
                
                if 'crc_stats' in block_info:
                    print(f"CRC32 checksums: {block_info['crc_stats']['valid']} valid, " +
                          f"{block_info['crc_stats']['invalid']} invalid, " +
                          f"{block_info['crc_stats']['missing']} missing")
                
                if block_info['block_size_bytes'] > 0:
                    print(f"Maximum block size: {block_info['block_size_bytes']} bytes")
                
                if block_info['actual_block_sizes']:
                    print("Actual block sizes:")
                    for block_idx, size in sorted(block_info['actual_block_sizes'].items()):
                        print(f"  Block {block_idx}: {size} bytes")
                
                # Update block size if it was read from the data
                if block_info['block_size_bytes'] > 0:
                    self.block_size_bytes = block_info['block_size_bytes']
                
                return True, block_info
            else:
                print("No blocks were successfully decoded")
                return False, {}
        except Exception as e:
            print(f"Error decoding DNA: {e}")
            return False, {}
    
    def reconstruct_file(self, block_files, output_file):
        """
        Reconstruct a file from multiple encoded blocks
        
        Args:
            block_files (list): List of paths to block CSV files
            output_file (str): Path to save the reconstructed file
            
        Returns:
            tuple: (Success status, File information)
        """
        # Dictionary to store blocks by their index
        blocks_by_index = {}
        file_info = {
            "total_blocks": 0,
            "total_bytes": 0,
            "pools_used": set()
        }
        
        # Process each block file
        for csv_file in block_files:
            # Determine Pool/Block from filename
            # Expected format: poolX_blockY.csv
            try:
                filename = os.path.basename(csv_file)
                pool_str = filename.split('_')[0].replace('pool', '')
                block_str = filename.split('_')[1].replace('block', '').split('.')[0]
                
                pool_idx = int(pool_str)
                block_idx = int(block_str)
            except:
                # If filename doesn't match expected pattern, try to extract from file content
                pool_idx = -1
                block_idx = -1
                
            # Set the current Pool/Block
            self.Pool = pool_idx
            self.Block = block_idx
            
            # Create temporary output file for this block
            temp_out_file = os.path.join(self.temp_dir, f"block_{block_idx}.bin")
            
            # Decode this block
            success, block_info = self.decode(csv_file, temp_out_file)
            
            if success:
                # Add to blocks dictionary
                file_info["pools_used"].add(pool_idx)
                
                if block_idx in blocks_by_index:
                    # We already have this block (from another pool)
                    # Keep the one with valid CRC if possible
                    existing_block = blocks_by_index[block_idx]
                    if ('crc_stats' in block_info and 
                        block_info['crc_stats']['valid'] > existing_block['info']['crc_stats']['valid']):
                        # This block has better CRC stats
                        blocks_by_index[block_idx] = {
                            "file": temp_out_file,
                            "info": block_info
                        }
                else:
                    # First time seeing this block
                    blocks_by_index[block_idx] = {
                        "file": temp_out_file,
                        "info": block_info
                    }
        
        # Update file info
        file_info["total_blocks"] = len(blocks_by_index)
        
        # Check if we have all blocks
        max_block_idx = max(blocks_by_index.keys()) if blocks_by_index else -1
        if max_block_idx >= 0 and all(i in blocks_by_index for i in range(max_block_idx + 1)):
            # We have all blocks - combine them
            with open(output_file, 'wb') as out_f:
                for block_idx in range(max_block_idx + 1):
                    block_data = blocks_by_index[block_idx]
                    with open(block_data["file"], 'rb') as in_f:
                        data = in_f.read()
                        out_f.write(data)
                        file_info["total_bytes"] += len(data)
            
            print(f"Successfully reconstructed file with {file_info['total_blocks']} blocks")
            print(f"Total size: {file_info['total_bytes']} bytes")
            print(f"Used {len(file_info['pools_used'])} pools: {sorted(file_info['pools_used'])}")
            
            return True, file_info
        else:
            print(f"Missing blocks. Have blocks: {sorted(blocks_by_index.keys())}")
            print(f"Expected blocks 0-{max_block_idx}")
            return False, file_info
    
    def get_all_pool_block_combinations(self):
        """
        Get all registered pool-block combinations
        
        Returns:
            list: List of (pool, block) tuples
        """
        return list(self.primer_registry.keys())
    
    def FilterSequence(self, sequence):
        """
        Filter the sequence to remove primers
        
        Args:
            sequence (str): DNA sequence
            
        Returns:
            str: Filtered sequence
        """
        # Get primers for the current Pool/Block
        forward_primer, reverse_primer = self.get_primers(self.Pool, self.Block)
        
        # Remove forward primer if present
        if forward_primer and sequence.startswith(forward_primer):
            sequence = sequence[len(forward_primer):]
        
        # Remove reverse primer if present
        if reverse_primer and sequence.endswith(reverse_primer):
            sequence = sequence[:-len(reverse_primer)]
        
        return sequence
    
    def ClassifySequence(self, sequence):
        """
        Determine which Pool and Block a sequence belongs to based on primers
        
        This method focuses only on Pool/Block assignment using primers,
        not on the placement of data within the file structure.
        
        Args:
            sequence (str): DNA sequence
            
        Returns:
            tuple: (Pool, Block, filtered_sequence)
        """
        oSeq = sequence
        cPool = -1
        cBlock = -1
        
        # Try all registered primer pairs
        for (pool, block), (forward_primer, reverse_primer) in self.primer_registry.items():
            if sequence.startswith(forward_primer):
                cPool = pool
                cBlock = block
                self.Pool = pool
                self.Block = block
                
                # Filter the sequence (remove primers)
                oSeq = sequence[len(forward_primer):]
                if oSeq.endswith(reverse_primer):
                    oSeq = oSeq[:-len(reverse_primer)]
                
                return cPool, cBlock, oSeq
        
        # If no registered primer matches, check if default primers are used
        if sequence.startswith(self.DEFAULT_FORWARD_PRIMER):
            # In this case, we can't determine Pool/Block from primers
            # but we can still return the filtered sequence
            oSeq = sequence[len(self.DEFAULT_FORWARD_PRIMER):]
            if oSeq.endswith(self.DEFAULT_REVERSE_PRIMER):
                oSeq = oSeq[:-len(self.DEFAULT_REVERSE_PRIMER)]
        
        # Return unknown Pool/Block with the filtered sequence
        # MolFS would need to use other mechanisms to determine Pool/Block
        return cPool, cBlock, oSeq
    
    def ClassifyByAddress(self, sequence):
        """
        Extract address information from a sequence for data placement
        
        This method focuses on getting the address information that determines 
        where the data payload should be placed in the file structure,
        independent of Pool/Block assignment.
        
        Args:
            sequence (str): DNA sequence (with or without primers)
            
        Returns:
            tuple: (address_int, address_bits, filtered_sequence)
        """
        # Check if there are known primers
        forward_primer, reverse_primer = self.get_primers(self.Pool, self.Block)
        
        # Extract core sequence without primers
        core_sequence = sequence
        if forward_primer and sequence.startswith(forward_primer):
            core_sequence = core_sequence[len(forward_primer):]
        if reverse_primer and core_sequence.endswith(reverse_primer):
            core_sequence = core_sequence[:-len(reverse_primer)]
            
        # Extract address from the first 19 nucleotides
        if len(core_sequence) >= 19:
            address_dna = core_sequence[:19]
            address_bits = ''.join(BASE_TO_BIT[base] for base in address_dna)
            address_int = int(address_bits, 2)
            
            return address_int, address_bits, core_sequence
        
        # If we can't find address information, return defaults
        return -1, "", core_sequence
    
    def ClassifySequence_SW(self, sequence):
        """
        Use Smith-Waterman algorithm to determine Pool and Block
        
        This is a placeholder for a future implementation using
        Smith-Waterman alignment for better error tolerance.
        
        Args:
            sequence (str): DNA sequence
            
        Returns:
            tuple: (Pool, Block, filtered_sequence)
        """
        # This is a placeholder for the Smith-Waterman implementation
        # For now, return unknown Pool/Block
        return -1, -1, sequence
    
    def getPrimerPair(self, cPool, cBlock):
        """
        Get the primer pair used for a specific Pool and Block
        
        Args:
            cPool (int): Pool number
            cBlock (int): Block number
            
        Returns:
            tuple: (forward_primer, reverse_primer)
        """
        return self.get_primers(cPool, cBlock)
    
    def exportSequences(self, in_file):
        """
        Read a block file and return the plain list of sequences
        
        Args:
            in_file (str): Path to the CSV file with DNA sequences
            
        Returns:
            list: List of DNA sequences
        """
        sequences = []
        
        # Read sequences from CSV
        with open(in_file, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile)
            header = next(reader)  # Skip header row
            
            # Find the column with DNA oligos
            oligo_idx = header.index("DNA Oligo") if "DNA Oligo" in header else 3
            
            for row in reader:
                sequences.append(row[oligo_idx])
        
        return sequences
