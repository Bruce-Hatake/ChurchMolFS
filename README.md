# ChurchMolFS
Updated old 2012 Church algorithm with CRC32 checksum, parameterized primers, file functionalities and integration with MolFS

# Church Algorithm for MolFS detailed README

This repository implements George Church's DNA storage algorithm for integration with the Molecular File System (MolFS). This implementation provides a complete solution for encoding binary data into DNA sequences, storing them in multiple pools with different primers, and decoding them back to the original files.

## Table of Contents

1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Module Descriptions](#module-descriptions)
    - [Encoder](#encoder-church_algorithm_encoderpy)
    - [Decoder](#decoder-church_algorithm_decoderpy)
    - [MolFS Interface](#molfs-interface-church_interfacepy)
    - [Test Script](#test-script-test_church_algorithmpy)
4. [DNA Oligo Structure](#dna-oligo-structure)
5. [Key Features](#key-features)
6. [Primer Configuration Guide](#primer-configuration-guide)
7. [Multi-File Management](#multi-file-management)
8. [Terminology](#terminology)
9. [Usage Examples](#usage-examples)
10. [Advanced Features](#advanced-features)
11. [Performance Considerations](#performance-considerations)
12. [Future Enhancements](#future-enhancements)

## Overview

The Church Algorithm, proposed by George Church et al. in "Next-Generation Digital Information Storage in DNA" (2012), offers a robust approach to storing digital information in DNA. This implementation adapts the algorithm for integration with the Molecular File System (MolFS), providing block-based storage, multi-pool distribution, and error detection.

## Architecture

The implementation consists of three main components with clear separation of concerns:

1. **Encoder**: Converts binary data to DNA sequences
2. **Decoder**: Converts DNA sequences back to binary data
3. **MolFS Interface**: Connects the encoder/decoder with MolFS, handling pool and block management

```
┌───────────────┐       ┌───────────────┐
│    MolFS      │◄─────►│   Interface   │
└───────────────┘       └───────┬───────┘
                              ▲ │
                              │ ▼
                     ┌────────┴────────┐
                     │                 │
               ┌─────┴─────┐     ┌─────┴─────┐
               │  Encoder  │     │  Decoder  │
               └───────────┘     └───────────┘
```

## Module Descriptions

### Encoder (`church_algorithm_encoder.py`)

The encoder module implements the Church DNA encoding scheme with several important features:

#### Key Functions:

- **encode_file_to_dna()**: Main entry point for encoding a file
- **bits_to_dna_with_crc32()**: Converts bits to DNA with CRC32 checksum
- **encode_bits_to_dna()**: Core function for bit-to-base conversion
- **avoid_homopolymer()**: Prevents runs of 4+ identical bases

#### Important Design Elements:

- **Bit-to-Base Mapping**: '0' is encoded as 'A' or 'C', '1' as 'G' or 'T'
- **Homopolymer Avoidance**: Automatically prevents runs of 4+ identical bases
- **Block-Based Processing**: Files are split into fixed-size blocks
- **CRC32 Checksums**: Each data payload has a 32-bit checksum for error detection

### Decoder (`church_algorithm_decoder.py`)

The decoder module implements the reverse process, converting DNA sequences back to the original binary data:

#### Key Functions:

- **decode_oligos_from_csv()**: Decodes DNA oligos from a CSV file
- **decode_oligos_from_fastq()**: Decodes DNA oligos from a FASTQ file
- **extract_oligo_parts()**: Extracts address, data, and checksum parts
- **verify_crc32_checksum()**: Validates data integrity

#### Important Design Elements:

- **Base-to-Bit Mapping**: 'A' and 'C' decode to '0', 'G' and 'T' decode to '1'
- **Error Detection**: CRC32 checksums are verified during decoding
- **Multi-Format Support**: Works with both CSV and FASTQ inputs
- **File Reconstruction**: Assembles decoded blocks in the correct order

### MolFS Interface (`church_interface.py`)

The interface module connects the encoder/decoder with MolFS, handling all pool and block management:

#### Key Methods:

- **encode_block()**: Encodes a specific file block with pool/block primers
- **encode_file()**: Encodes an entire file using a distribution strategy
- **decode()**: Decodes DNA sequences back to the original file
- **reconstruct_file()**: Reconstructs a file from multiple blocks/pools
- **ClassifySequence()**: Identifies pool/block from primers

#### Important Design Elements:

- **Primer Registry**: Stores custom primers for each pool/block combination
- **Pool/Block Management**: Maps blocks to pools with flexible strategies
- **Redundancy Support**: Same block can be stored in multiple pools
- **MolFS Compatibility**: Implements the required MolFS interface methods

### Test Script (`test_church_algorithm.py`)

The test script verifies the implementation with real-world scenarios:

- **Basic Functionality**: Simple encode/decode flow
- **Multi-Block**: File splitting across blocks
- **Redundant Storage**: Block replication across pools with recovery

## DNA Oligo Structure

Each DNA oligo has the following structure:

```
┌─────────────────┬─────────────┬────────────────┬───────────────┬─────────────────┐
│ Forward Primer  │   Address   │  Data Payload  │ CRC32 Checksum│ Reverse Primer  │
│    (22 nt)      │   (19 nt)   │    (96 nt)     │   (32 nt)     │    (22 nt)      │
└─────────────────┴─────────────┴────────────────┴───────────────┴─────────────────┘
```

1. **Forward Primer** (22 nt): Used for PCR and pool/block identification
2. **Address** (19 nt): Uniquely identifies the position within the block
3. **Data Payload** (96 nt): Carries 96 bits (12 bytes) of actual data
4. **CRC32 Checksum** (32 nt): Provides error detection for the data payload
5. **Reverse Primer** (22 nt): Used for PCR and reverse complementing

## Key Features

### 1. Block-Based Storage

The implementation divides files into fixed-size blocks (default 5KB), which provides several advantages:

- **Scalability**: Handles files of any size
- **Distribution**: Flexible assignment of blocks to pools
- **Redundancy**: Can store the same block in multiple pools
- **Recovery**: Can reconstruct files even if some blocks are missing

Example:
```python
molfs_dev = MolFSDev()
molfs_dev.set_block_size(5 * 1024)  # Set block size to 5KB
```

### 2. Redundant Storage

The same block can be stored in multiple pools with different primers, providing data durability:

```python
# Store Block 1 in both Pool 1 and Pool 2
molfs_dev.register_primers(1, 1, "CTACACGACACTTTTCCGATCT", "AGATCGGAAGAGCGGAACAGCA")
molfs_dev.register_primers(2, 1, "CTACACGACGCTTAACCGATCT", "AGATCGGAAGAGCGGATCAGCA")

# Define distribution strategy
def redundant_strategy(block_idx, total_blocks):
    if block_idx == 1:
        return [1, 2]  # Block 1 goes to both pools
    else:
        return [1]     # Other blocks just go to Pool 1
```

### 3. Custom Distribution Strategies

The implementation allows custom distribution strategies for assigning blocks to pools:

```python
# Example: Store even blocks in Pool 1, odd blocks in Pool 2
def even_odd_strategy(block_idx, total_blocks):
    if block_idx % 2 == 0:
        return [1]  # Even blocks to Pool 1
    else:
        return [2]  # Odd blocks to Pool 2
```

### 4. Error Detection

CRC32 checksums provide error detection for each data payload:

- Checksums are calculated during encoding
- Verified during decoding
- Invalid data is automatically discarded
- Statistics are provided on checksum validation

### 5. Multi-Format Support

The implementation supports both CSV and FASTQ input formats:

- **CSV**: For processed data with metadata
- **FASTQ**: For raw sequencing data

## Primer Configuration Guide

Primers play a crucial role in this implementation as they determine the Pool/Block assignment and enable PCR amplification. This section provides detailed guidance on configuring primers for different scenarios.

### Important Requirements

1. **Primer Length**: All primers must be **exactly 22 bases long**
2. **Primer Diversity**: Primers should be sufficiently different to avoid cross-hybridization
3. **Balanced GC Content**: Aim for 40-60% GC content for optimal PCR performance

### Configuring Multiple Blocks in the Same Pool

When a pool contains multiple different blocks, each block should have a unique primer pair:

```python
# Pool 1 with 3 different blocks
molfs_dev.register_primers(1, 0, "CTACACGACGCTCTTCCGATCT", "AGATCGGAAGAGCGGTTCAGCA")  # Block 0
molfs_dev.register_primers(1, 1, "CTACACGACACTTTTCCGATCT", "AGATCGGAAGAGCGGAACAGCA")  # Block 1
molfs_dev.register_primers(1, 2, "CTACACGACGCTAATCCGATCT", "AGATCGGAAGAGCGGGTCAGCA")  # Block 2
```

This configuration allows the system to:
- Distinguish between blocks within the same pool
- Selectively amplify specific blocks during PCR
- Correctly classify sequences during decoding

### Configuring the Same Block Across Different Pools

For redundancy, the same block can be stored in multiple pools, each with a different primer pair:

```python
# Block 1 in both Pool 1 and Pool 2 with different primers
molfs_dev.register_primers(1, 1, "CTACACGACACTTTTCCGATCT", "AGATCGGAAGAGCGGAACAGCA")  # Pool 1, Block 1
molfs_dev.register_primers(2, 1, "CTACACGACGCTTAACCGATCT", "AGATCGGAAGAGCGGATCAGCA")  # Pool 2, Block 1
```

This configuration provides:
- Redundant storage of the same block
- Ability to recover if one pool is damaged
- Distinct identification of which pool a sequence belongs to

### Designing Diverse Primers

When designing primers for many Pool/Block combinations, follow these guidelines:

1. **Maintain Fixed Ends**: Keep the ends consistent for PCR compatibility
   ```
   CTACACGAC[variable 6nt]CCGATCT
   ```
   This ensures the total primer length is always exactly 22 nucleotides.

2. **Vary the Middle Section**: Create diversity in the central 6-nucleotide portion
   ```
   CTACACGACGCTCTTCCGATCT  # Original
   CTACACGACACTTTTCCGATCT  # Variation 1
   CTACACGACGCTTAACCGATCT  # Variation 2
   CTACACGACGCTAATCCGATCT  # Variation 3
   ```

3. **Ensure Sufficient Hamming Distance**: At least 3-4 base differences between primers

4. **Consider Professional Primer Design Tools**: For production use, consider specialized primer design software that accounts for:
   - Melting temperature (Tm)
   - Secondary structure formation
   - Primer-dimer avoidance
   - GC content optimization

5. **Example of a Simple Systematic Approach**: 
   This is just one possible method, not necessarily optimal for all cases:
   ```python
   # An example approach to generate a primer for a specific pool and block
   def generate_primer(pool, block):
       # Create unique variations based on pool and block IDs
       var1 = "ACGT"[pool % 4]
       var2 = "ACGT"[block % 4]
       var3 = "ACGT"[(pool + block) % 4]
       
       # Create a unique 6nt middle section
       middle = f"G{var1}{var2}{var3}C"  # Results in 6nt
       
       # Combine with standard ends (total: 22nt)
       return f"CTACACGAC{middle}CCGATCT"  # 9 + 6 + 7 = 22nt
   ```

## Multi-File Management

The current implementation focuses on the storage and retrieval of DNA-encoded data blocks, but it doesn't inherently include a file identification system within the DNA oligos themselves. This section explains how to handle multiple files in the MolFS system.

### Address Scope and Limitations

The 19-bit address in each DNA oligo identifies:
- The position within a specific block (which oligo it is)
- It does NOT identify which file the block belongs to

This means that when storing multiple files, additional tracking mechanisms are needed.

### File-to-Block Mapping

For a complete file system, MolFS should implement a file-to-block mapping system:

1. **File Registry**: MolFS should maintain a database/registry that tracks:
   - Which files are stored in the system
   - Which blocks belong to each file
   - Which pools contain each block

2. **Block Assignment**: When encoding a file, MolFS should:
   - Assign unique block IDs to the file
   - Register these mappings in its internal database
   - Use the appropriate primers for each block/pool combination

3. **File Reconstruction**: When retrieving a file, MolFS should:
   - Look up which blocks constitute the file
   - Retrieve those blocks from their respective pools
   - Use the `reconstruct_file` method to combine the blocks

### Example Implementation

```python
# MolFS file-to-block mapping database (simplified example)
file_registry = {
    "document1.pdf": {
        "blocks": [0, 1],
        "pools": {0: [1], 1: [1, 2]}  # Block 0 in Pool 1, Block 1 in Pools 1 and 2
    },
    "image.jpg": {
        "blocks": [2, 3],
        "pools": {2: [2], 3: [2]}  # Blocks 2 and 3 both in Pool 2
    }
}

# To retrieve document1.pdf
file_info = file_registry["document1.pdf"]
block_files = []

# Collect all block files (using the first available pool for each block)
for block_id in file_info["blocks"]:
    pool_id = file_info["pools"][block_id][0]  # Use first available pool
    block_file = f"pool{pool_id}_block{block_id}.csv"
    block_files.append(block_file)

# Reconstruct the file
molfs_dev.reconstruct_file(block_files, "document1.pdf")
```

### Alternative Approaches

For more advanced implementations, consider:

1. **Metadata Blocks**: Each file could have a special metadata block that:
   - Lists all blocks belonging to the file
   - Contains file properties (name, size, creation date, etc.)
   - Has a higher replication factor for durability

2. **Embedded File ID**: Future versions could modify the DNA structure to include:
   - File identifier in the address portion
   - Version information for files that change over time
   - Content-based addressing for deduplication

This separation of concerns:
- Keeps the DNA oligo structure simple and efficient
- Allows MolFS to implement sophisticated file management
- Enables features like versioning, deduplication, and access control

## Terminology

To avoid confusion, it's important to understand the terminology used:

- **MolFS Block**: A larger organizational unit (e.g., 5KB) that MolFS uses to distribute data across pools
- **DNA Oligo**: An individual DNA sequence containing address, data payload, and checksum
- **Pool**: A collection of blocks with unique primers
- **Primer**: A DNA sequence used for PCR and pool/block identification

## Usage Examples

### Basic Usage

```python
from church_interface import MolFSDev

# Initialize
molfs_dev = MolFSDev()
molfs_dev.set_block_size(5 * 1024)  # 5KB blocks

# Register primers (must be exactly 22 bases)
molfs_dev.register_primers(1, 0, "CTACACGACGCTCTTCCGATCT", "AGATCGGAAGAGCGGTTCAGCA")

# Set current pool/block
molfs_dev.Pool = 1
molfs_dev.Block = 0

# Encode a file
molfs_dev.encode("input.bin", "output.csv")

# Decode a file
molfs_dev.decode("output.csv", "reconstructed.bin")
```

### Multi-Pool Distribution

```python
# Define a distribution strategy
def my_strategy(block_idx, total_blocks):
    if block_idx < total_blocks // 2:
        return [1]  # First half of blocks to Pool 1
    else:
        return [2]  # Second half of blocks to Pool 2

# Register primers for all pools and blocks
molfs_dev.register_primers(1, 0, "CTACACGACGCTCTTCCGATCT", "AGATCGGAAGAGCGGTTCAGCA")
molfs_dev.register_primers(1, 1, "CTACACGACACTTTTCCGATCT", "AGATCGGAAGAGCGGAACAGCA")
molfs_dev.register_primers(2, 2, "CTACACGACGCTAATCCGATCT", "AGATCGGAAGAGCGGGTCAGCA")
molfs_dev.register_primers(2, 3, "CTACACGACGCTTAACCGATCT", "AGATCGGAAGAGCGGATCAGCA")

# Encode file across multiple pools
block_distribution = molfs_dev.encode_file("input.bin", "output_dir", my_strategy)

# Get all block files
block_files = [info["file"] for _, info in block_distribution.items()]

# Reconstruct the file
molfs_dev.reconstruct_file(block_files, "reconstructed.bin")
```

## Advanced Features

### 1. Pool Failure Recovery

The system can recover from pool failures by using redundant blocks:

```python
# Get all files except those from the failed pool
available_files = [info["file"] for (pool, _), info in block_distribution.items() 
                  if pool != failed_pool]

# Reconstruct with available blocks
molfs_dev.reconstruct_file(available_files, "reconstructed.bin")
```

### 2. CRC32 Statistics

The decoder provides detailed statistics on CRC32 checksum validation:

```python
_, block_info = molfs_dev.decode("encoded.csv", "decoded.bin")
print(f"CRC32 checksums: {block_info['crc_stats']['valid']} valid, " +
      f"{block_info['crc_stats']['invalid']} invalid, " +
      f"{block_info['crc_stats']['missing']} missing")
```

### 3. Block Size Reporting

The system tracks actual block sizes, which is especially important for the last block that might be smaller:

```python
for block_idx, size in sorted(block_info['actual_block_sizes'].items()):
    print(f"Block {block_idx}: {size} bytes")
```

## Performance Considerations

1. **DNA Length**: Each oligo can store 96 bits (12 bytes) of payload data
2. **Oligo Count**: A 5KB block requires ~450 oligos, each with its own address
3. **CRC32 Overhead**: The CRC32 checksum adds 32 bits of overhead per oligo
4. **Primer Length**: Primers add 44nt (22nt each) to each oligo

## Future Enhancements

This implementation provides a solid foundation that could be extended with:

1. **Consensus Building**: For handling sequencing errors in real-world applications
2. **Error Correction Codes**: Reed-Solomon or LDPC codes for error correction
3. **Compression**: Pre-encoding compression to increase storage density
4. **Encryption**: Pre-encoding encryption for data security
5. **Integrated File ID**: Embedding file identifiers in the DNA structure

## Integration Notes for MolFS

When integrating with MolFS, consider the following:

1. **Primer Design**: Design primers that are sufficiently different for reliable discrimination
2. **Block Size**: Choose block sizes based on your storage strategy and file characteristics
3. **Redundancy Strategy**: Consider criticality of data when deciding on redundancy level
4. **File Tracking**: Implement a robust file-to-block mapping system
5. **Addressing**: The 19-bit address allows for up to 524,288 oligos per block

For optimal integration, you might want to:

1. Create a higher-level module that manages the distribution of files across pools
2. Develop a database to track which blocks are stored in which pools
3. Implement recovery strategies based on the redundancy level
4. Consider adding a consensus building module for handling sequencing errors

---

Feel free to reach out if you need any clarification or have questions about the implementation.
