%%writefile {save_dir}church_algorithm_decoder.py


import csv
import zlib
import os

# Mapping bases back to bits
BASE_TO_BIT = {'A': '0', 'C': '0', 'G': '1', 'T': '1'}

# Default primer sequences
DEFAULT_FORWARD_PRIMER = "CTACACGACGCTCTTCCGATCT"
DEFAULT_REVERSE_PRIMER = "AGATCGGAAGAGCGGTTCAGCA"

def decode_dna_to_bits(dna_sequence):
    """
    Decode a DNA sequence back into a bitstring.
    
    Args:
        dna_sequence (str): Input DNA sequence.
        
    Returns:
        str: Decoded bitstring.
    """
    return ''.join(BASE_TO_BIT[base] for base in dna_sequence)

def verify_crc32_checksum(data_bits, checksum_bits):
    """
    Verify that the CRC32 checksum of the data matches the provided checksum.
    
    Args:
        data_bits (str): Data bits to verify.
        checksum_bits (str): Expected CRC32 checksum as bits.
        
    Returns:
        bool: True if checksum matches, False otherwise.
    """
    # Convert data bits to bytes
    byte_data = int(data_bits, 2).to_bytes((len(data_bits) + 7) // 8, byteorder='big')
    
    # Calculate checksum
    calculated_checksum = zlib.crc32(byte_data)
    calculated_checksum_bits = f"{calculated_checksum:032b}"
    
    return calculated_checksum_bits == checksum_bits

def extract_oligo_parts(oligo, forward_primer=None, reverse_primer=None):
    """
    Extract address, data, and CRC32 parts from an oligo.
    
    Args:
        oligo (str): Full DNA oligo sequence.
        forward_primer (str, optional): Forward primer sequence.
        reverse_primer (str, optional): Reverse primer sequence.
        
    Returns:
        tuple: (address_dna, data_dna, crc32_dna)
    """
    # If primers are provided, remove them
    core_sequence = oligo
    
    if forward_primer and oligo.startswith(forward_primer):
        core_sequence = core_sequence[len(forward_primer):]
    
    if reverse_primer and core_sequence.endswith(reverse_primer):
        core_sequence = core_sequence[:-len(reverse_primer)]
    
    # Extract parts based on known positions
    # Address: 19 nt, Data: 96 nt, CRC32: rest
    if len(core_sequence) >= 19 + 96:
        address_dna = core_sequence[:19]
        data_dna = core_sequence[19:19+96]
        crc32_dna = core_sequence[19+96:]
        return address_dna, data_dna, crc32_dna
    
    # If the sequence is too short, return None
    return None, None, None

def decode_oligos_from_csv(csv_file_path, forward_primer=None, reverse_primer=None, output_file=None):
    """
    Decode DNA oligos from a CSV file back into the original data.
    
    Args:
        csv_file_path (str): Path to the CSV file containing DNA oligos.
        forward_primer (str, optional): Forward primer sequence.
        reverse_primer (str, optional): Reverse primer sequence.
        output_file (str, optional): Path to save the decoded binary file.
        
    Returns:
        tuple: (Decoded blocks dict, Block information dict)
            - decoded_blocks: Dictionary mapping (block_index, address) to decoded data bits
            - block_info: Dictionary with block size and count information
    """
    decoded_blocks = {}
    errors = []
    block_info = {
        "block_size_bytes": 0,      # Maximum block size
        "block_count": 0,           # Number of blocks
        "total_file_size": 0,       # Total size of the original file
        "actual_block_sizes": {},    # Actual size of each block (key: block_index)
        "crc_stats": {"valid": 0, "invalid": 0, "missing": 0}  # CRC checksum statistics
    }
    
    # Read oligos from CSV
    with open(csv_file_path, 'r', newline='') as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader)  # Skip header row
        
        # Determine which columns contain the relevant data
        oligo_idx = header.index("DNA Oligo") if "DNA Oligo" in header else -1
        block_idx_col = header.index("Block Index") if "Block Index" in header else -1
        block_size_col = header.index("Block Size (Bytes)") if "Block Size (Bytes)" in header else -1
        actual_block_size_col = header.index("Actual Block Size (Bytes)") if "Actual Block Size (Bytes)" in header else -1
        total_file_size_col = header.index("Total File Size (Bytes)") if "Total File Size (Bytes)" in header else -1
        
        if oligo_idx == -1:
            raise ValueError("CSV file does not contain a 'DNA Oligo' column")
        
        block_indices = set()
        
        for row_num, row in enumerate(reader, 2):  # Start at row 2 (accounting for header)
            try:
                oligo = row[oligo_idx].strip()
                
                # Get block index if available
                block_index = int(row[block_idx_col]) if block_idx_col != -1 and row[block_idx_col] else 0
                block_indices.add(block_index)
                
                # Get block size if available
                if block_size_col != -1 and row[block_size_col]:
                    block_info["block_size_bytes"] = int(row[block_size_col])
                
                # Get actual block size if available
                if actual_block_size_col != -1 and row[actual_block_size_col]:
                    actual_block_size = int(row[actual_block_size_col])
                    block_info["actual_block_sizes"][block_index] = actual_block_size
                
                # Get total file size if available
                if total_file_size_col != -1 and row[total_file_size_col]:
                    block_info["total_file_size"] = int(row[total_file_size_col])
                
                # Extract parts from the oligo
                address_dna, data_dna, crc32_dna = extract_oligo_parts(oligo, forward_primer, reverse_primer)
                
                if address_dna is None or data_dna is None:
                    errors.append(f"Row {row_num}: Failed to extract parts from oligo")
                    continue
                
                # Decode DNA to bits
                address_bits = decode_dna_to_bits(address_dna)
                data_bits = decode_dna_to_bits(data_dna)
                
                # Handle CRC32 checksum if available
                crc_valid = False
                if crc32_dna and len(crc32_dna) >= 32:
                    crc32_bits = decode_dna_to_bits(crc32_dna[:32])
                    crc_valid = verify_crc32_checksum(data_bits, crc32_bits)
                    if crc_valid:
                        block_info["crc_stats"]["valid"] += 1
                    else:
                        block_info["crc_stats"]["invalid"] += 1
                        errors.append(f"Row {row_num}: CRC32 checksum verification failed")
                        continue  # Skip this oligo if CRC is invalid
                else:
                    block_info["crc_stats"]["missing"] += 1
                
                # Convert address to integer for sorting
                address = int(address_bits, 2)
                
                # Store with block information
                decoded_blocks[(block_index, address)] = data_bits
                
            except Exception as e:
                errors.append(f"Row {row_num}: Error processing - {str(e)}")
        
        # Update block count
        block_info["block_count"] = len(block_indices)
    
    # If there were errors, report them
    if errors:
        print(f"Encountered {len(errors)} errors during decoding:")
        for error in errors[:10]:  # Show first 10 errors
            print(f"  - {error}")
        if len(errors) > 10:
            print(f"  ... and {len(errors) - 10} more errors")
    
    # Write decoded data to file if specified
    if output_file and decoded_blocks:
        # Sort by block index and then by address within each block
        sorted_keys = sorted(decoded_blocks.keys())
        
        # Concatenate all data bits in the correct order
        reconstructed_bitstream = ""
        for key in sorted_keys:
            reconstructed_bitstream += decoded_blocks[key]
        
        # If we know the total file size, trim the reconstructed bitstream
        if block_info["total_file_size"] > 0:
            # Convert file size to bits
            total_file_bits = block_info["total_file_size"] * 8
            reconstructed_bitstream = reconstructed_bitstream[:total_file_bits]
        
        # Convert bitstream to bytes
        bit_length = len(reconstructed_bitstream)
        byte_length = bit_length // 8
        
        # Only process complete bytes
        if byte_length > 0:
            byte_data = int(reconstructed_bitstream[:byte_length*8], 2).to_bytes(byte_length, byteorder='big')
            
            with open(output_file, 'wb') as f:
                f.write(byte_data)
            
            print(f"Decoded data written to {output_file} ({byte_length} bytes)")
    
    return decoded_blocks, block_info

def decode_oligos_from_fastq(fastq_file, forward_primer=None, reverse_primer=None, output_file=None):
    """
    Basic decoder for DNA oligos from a FASTQ file - no error correction.
    This is a simplified version that processes each read independently.
    For a real-world implementation, you would want to use a separate module
    for consensus building and error correction.
    
    Args:
        fastq_file (str): Path to the FASTQ file containing DNA sequences.
        forward_primer (str, optional): Forward primer sequence.
        reverse_primer (str, optional): Reverse primer sequence.
        output_file (str, optional): Path to save the decoded binary file.
        
    Returns:
        tuple: (Decoded blocks dict, Block information dict)
    """
    # Initialize data structures
    decoded_blocks = {}
    block_info = {
        "block_size_bytes": 0,
        "block_count": 0,
        "total_file_size": 0,
        "actual_block_sizes": {},
        "crc_stats": {"valid": 0, "invalid": 0, "missing": 0}
    }
    
    # Read sequences from FASTQ file
    sequences = []
    with open(fastq_file, 'r') as f:
        line_count = 0
        for line in f:
            line_count += 1
            if line_count % 4 == 2:  # FASTQ format: every 4th line is sequence
                sequences.append(line.strip())
    
    print(f"Read {len(sequences)} sequences from {fastq_file}")
    
    # Process each sequence independently
    for seq in sequences:
        try:
            # Extract parts from the sequence
            address_dna, data_dna, crc32_dna = extract_oligo_parts(seq, forward_primer, reverse_primer)
            
            if address_dna is None or data_dna is None:
                continue
            
            # Decode DNA to bits
            address_bits = decode_dna_to_bits(address_dna)
            data_bits = decode_dna_to_bits(data_dna)
            
            # Check CRC32 if available
            crc_valid = False
            if crc32_dna and len(crc32_dna) >= 32:
                crc32_bits = decode_dna_to_bits(crc32_dna[:32])
                crc_valid = verify_crc32_checksum(data_bits, crc32_bits)
                if crc_valid:
                    block_info["crc_stats"]["valid"] += 1
                else:
                    block_info["crc_stats"]["invalid"] += 1
                    continue  # Skip if CRC is invalid
            else:
                block_info["crc_stats"]["missing"] += 1
            
            # Convert address to integer and use default block 0
            address = int(address_bits, 2)
            block_index = 0  # Default block
            
            # Store data bits
            key = (block_index, address)
            decoded_blocks[key] = data_bits
            
        except Exception:
            # Skip sequences that can't be processed
            continue
    
    # Update block count and file size
    block_indices = set(block for block, _ in decoded_blocks.keys())
    block_info["block_count"] = len(block_indices)
    
    if decoded_blocks:
        # Estimate total file size from decoded blocks
        block_info["total_file_size"] = sum(len(bits) for bits in decoded_blocks.values()) // 8
    
    # Write decoded data to file if specified
    if output_file and decoded_blocks:
        # Sort by block index and then by address within each block
        sorted_keys = sorted(decoded_blocks.keys())
        
        # Concatenate all data bits in the correct order
        reconstructed_bitstream = ""
        for key in sorted_keys:
            reconstructed_bitstream += decoded_blocks[key]
        
        # Convert bitstream to bytes
        bit_length = len(reconstructed_bitstream)
        byte_length = bit_length // 8
        
        # Only process complete bytes
        if byte_length > 0:
            byte_data = int(reconstructed_bitstream[:byte_length*8], 2).to_bytes(byte_length, byteorder='big')
            
            with open(output_file, 'wb') as f:
                f.write(byte_data)
            
            print(f"Decoded data written to {output_file} ({byte_length} bytes)")
    
    return decoded_blocks, block_info

def bitstream_to_bytes(bitstream):
    """
    Convert a bitstream to bytes.
    
    Args:
        bitstream (str): Bitstream as a string of '0' and '1'.
        
    Returns:
        bytes: Converted bytes.
    """
    # Make sure the bitstream length is a multiple of 8
    bit_length = len(bitstream)
    byte_length = bit_length // 8
    
    if byte_length == 0:
        return b''
    
    # Convert bits to bytes
    return int(bitstream[:byte_length*8], 2).to_bytes(byte_length, byteorder='big')
