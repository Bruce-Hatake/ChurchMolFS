%%writefile {save_dir}church_algorithm_encoder.py


import csv
import zlib  # For CRC32 checksum
import os    # For checking file size
import random

# Mapping for bits to bases (with synonyms)
BIT_TO_BASE = {'0': ['A', 'C'], '1': ['G', 'T']}

# Universal primers (default values)
DEFAULT_FORWARD_PRIMER = "CTACACGACGCTCTTCCGATCT"
DEFAULT_REVERSE_PRIMER = "AGATCGGAAGAGCGGTTCAGCA"

# Avoid homopolymer runs of 4 bases
HOMOPOLYMER_LIMIT = 3

def avoid_homopolymer(last_bases, base):
    """Ensure no homopolymer runs of 4 bases."""
    if len(last_bases) >= HOMOPOLYMER_LIMIT and all(b == base for b in last_bases[-HOMOPOLYMER_LIMIT:]):
        synonyms = {'A': 'C', 'C': 'A', 'G': 'T', 'T': 'G'}
        return synonyms[base]
    return base

def encode_bits_to_dna(bits, last_bases):
    """Convert a sequence of bits into a DNA sequence."""
    dna_seq = ""
    for bit in bits:
        base = random.choice(BIT_TO_BASE[bit])  # Randomly choose between synonyms
        base = avoid_homopolymer(last_bases, base)  # Avoid homopolymer runs
        last_bases.append(base)
        if len(last_bases) > HOMOPOLYMER_LIMIT:
            last_bases.pop(0)
        dna_seq += base
    return dna_seq

def calculate_crc32_checksum(data_bits):
    """Calculate a 32-bit CRC checksum for the given data bits."""
    # Convert bits to bytes
    byte_data = int(data_bits, 2).to_bytes((len(data_bits) + 7) // 8, byteorder='big')
    checksum = zlib.crc32(byte_data)  # 32-bit CRC checksum
    return f"{checksum:032b}"  # Return checksum as a 32-bit binary string

def bits_to_dna_with_crc32(bits, forward_primer, reverse_primer, segment_length=12, max_address=524288):
    """
    Encode a bitstream string into DNA oligos with CRC32 checksum.

    Args:
        bits (str): Input bitstream.
        forward_primer (str): Forward primer sequence.
        reverse_primer (str): Reverse primer sequence.
        segment_length (int): Length of data in bytes per segment.
        max_address (int): Maximum number of addresses.

    Returns:
        list: List of tuples (address, data_bits, crc32_checksum, DNA oligo).
    """
    oligos_with_details = []
    address = 1  # Start address at 1 (0000000000000000001)
    last_bases = []

    # Validate input
    if not all(bit in "01" for bit in bits):
        raise ValueError("Input contains invalid characters. Only '0' and '1' are allowed.")

    # Process bitstream into DNA oligos
    for i in range(0, len(bits), segment_length * 8):  # Segment length is in bytes (8 bits/byte)
        if address > max_address:
            break  # Stop if max address is exceeded

        # Encode the 19-bit address
        address_bits = f"{address:019b}"  # Convert address to 19-bit binary
        address_dna = encode_bits_to_dna(address_bits, last_bases.copy())  # Use a copy to avoid affecting next encodings

        # Get the 96-bit data block
        data_bits = bits[i:i + segment_length * 8].ljust(96, '0')  # Ensure 96 bits with padding if needed

        # Calculate CRC32 checksum for the data block
        crc32_checksum_bits = calculate_crc32_checksum(data_bits)

        # Encode data into DNA
        data_dna = encode_bits_to_dna(data_bits, last_bases.copy())
        
        # Encode CRC32 into DNA
        crc32_dna = encode_bits_to_dna(crc32_checksum_bits, last_bases.copy())

        # Create the core oligo (address + data + CRC32)
        core_oligo = address_dna + data_dna + crc32_dna
        
        # Add primers if provided
        complete_oligo = core_oligo
        if forward_primer:
            complete_oligo = forward_primer + complete_oligo
        if reverse_primer:
            complete_oligo = complete_oligo + reverse_primer
            
        oligos_with_details.append((address_bits, data_bits, crc32_checksum_bits, complete_oligo))

        address += 1  # Increment address for the next segment

    return oligos_with_details

def read_binary_file(file_path):
    """
    Read a binary file and convert its contents to a bitstream.

    Args:
        file_path (str): Path to the binary file.

    Returns:
        tuple: (Bitstream as a string of '0' and '1', file size in bytes).
    """
    with open(file_path, 'rb') as file:
        binary_data = file.read()
    file_size = len(binary_data)
    return ''.join(f"{byte:08b}" for byte in binary_data), file_size  # Convert bytes to a string of bits

def save_to_csv_with_blocks(file_path, oligos_with_details, total_data_size, block_size_bytes):
    """
    Save DNA oligos and metadata to a CSV file with block information.
    
    Args:
        file_path (str): Path to the output CSV file.
        oligos_with_details (list): List of tuples (block_index, address_bits, data_bits, crc32_checksum, DNA oligo, actual_block_bytes).
        total_data_size (int): Total input data size in bytes.
        block_size_bytes (int): Maximum size of each block in bytes.
    """
    with open(file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "Block Index", 
            "Address", 
            "Data Bits", 
            "CRC32 Checksum", 
            "DNA Oligo", 
            "Block Size (Bytes)", 
            "Actual Block Size (Bytes)",
            "Total File Size (Bytes)"
        ])
        
        for block_index, address_bits, data_bits, crc32_checksum, oligo, actual_block_bytes in oligos_with_details:
            writer.writerow([
                block_index, 
                address_bits, 
                data_bits, 
                crc32_checksum, 
                oligo, 
                block_size_bytes,
                actual_block_bytes,
                total_data_size
            ])

def encode_file_to_dna(input_file_path, output_csv_path, forward_primer=None, reverse_primer=None, 
                   check_file_size=False, max_file_size=100*1024, block_size_bytes=1024):
    """
    Main function to encode a binary file into DNA sequences and save to CSV.
    
    Args:
        input_file_path (str): Path to the input binary file.
        output_csv_path (str): Path to save the output CSV file.
        forward_primer (str, optional): Forward primer sequence (22 bases).
        reverse_primer (str, optional): Reverse primer sequence (22 bases).
        check_file_size (bool, optional): Whether to check if file size exceeds maximum.
        max_file_size (int, optional): Maximum allowed file size in bytes.
        block_size_bytes (int, optional): Size of each DNA data block in bytes.
        
    Returns:
        tuple: (Number of oligos generated, Input file size in bytes, Number of blocks)
    """
    # Use default primers if not provided
    if forward_primer is None:
        forward_primer = DEFAULT_FORWARD_PRIMER
    if reverse_primer is None:
        reverse_primer = DEFAULT_REVERSE_PRIMER
        
    # Ensure primers are exactly 22 bases
    if forward_primer and len(forward_primer) != 22:
        raise ValueError("Forward primer must be exactly 22 bases long.")
    if reverse_primer and len(reverse_primer) != 22:
        raise ValueError("Reverse primer must be exactly 22 bases long.")
    
    # Check the file size if needed
    file_size = os.path.getsize(input_file_path)
    if check_file_size and file_size > max_file_size:
        raise ValueError(f"File size exceeds the maximum allowed size of {max_file_size / 1024} KB.")
    
    # Read the binary file and convert to bitstream
    bits, total_data_size = read_binary_file(input_file_path)
    
    # Calculate parameters for fixed block sizes
    # Each DNA oligo can store 12 bytes (96 bits) of data
    bytes_per_oligo = 12  # 96 bits / 8 bits per byte
    
    # Calculate how many oligos we need per block
    oligos_per_block = (block_size_bytes + bytes_per_oligo - 1) // bytes_per_oligo
    
    # Calculate total number of blocks needed
    total_bytes = len(bits) // 8
    num_blocks = (total_bytes + block_size_bytes - 1) // block_size_bytes
    
    # Create a list to store all oligos for all blocks
    all_oligos_with_details = []
    
    # Process the bitstream in block-sized chunks
    for block_index in range(num_blocks):
        # Calculate start and end positions for this block
        start_bit = block_index * block_size_bytes * 8
        end_bit = min(start_bit + block_size_bytes * 8, len(bits))
        block_bits = bits[start_bit:end_bit]
        
        # Calculate actual bytes in this block
        actual_block_bytes = (len(block_bits) + 7) // 8  # Round up to nearest byte
        
        # For partial blocks (typically the last block), calculate how many oligos we actually need
        actual_oligos_needed = (actual_block_bytes + bytes_per_oligo - 1) // bytes_per_oligo
        max_address_for_block = min(oligos_per_block, actual_oligos_needed)
        
        # Encode this block into DNA oligos
        block_oligos = bits_to_dna_with_crc32(
            block_bits, 
            forward_primer, 
            reverse_primer,
            segment_length=bytes_per_oligo,
            max_address=max_address_for_block
        )
        
        # Add block information to each oligo
        for i, (address_bits, data_bits, crc32_bits, oligo) in enumerate(block_oligos):
            # Add block index and actual block size to the tuple for tracking
            all_oligos_with_details.append((
                block_index, 
                address_bits, 
                data_bits, 
                crc32_bits, 
                oligo, 
                actual_block_bytes  # Store the actual size of this block
            ))
    
    # Save to a CSV file with block information
    save_to_csv_with_blocks(output_csv_path, all_oligos_with_details, total_data_size, block_size_bytes)
    
    return len(all_oligos_with_details), total_data_size, num_blocks
