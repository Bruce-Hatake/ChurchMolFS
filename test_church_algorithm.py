%%writefile {save_dir}test_church_algorithm.py


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script for the Church Algorithm implementation with MolFS integration

This script performs a simple test of the encoding and decoding functionality
to verify that the implementation works correctly.
"""

import os
import sys
import shutil
import random
import hashlib
from church_interface import MolFSDev

# Set base path for test files - adjust as needed
BASE_PATH = "/mnt/c/Users/ParkJ/Downloads/church_test/"

def create_test_file(path, size):
    """Create a test file with recognizable data pattern."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'wb') as f:
        # Create a pattern: A block of 'A', a block of 'B', etc.
        block_size = 256
        for i in range(0, size, block_size):
            char = bytes([65 + (i // block_size) % 26])  # A, B, C, ...
            f.write(char * min(block_size, size - i))

def verify_file_integrity(original_file, reconstructed_file):
    """Verify that the reconstructed file matches the original."""
    if not os.path.exists(reconstructed_file):
        print(f"Error: Reconstructed file {reconstructed_file} does not exist.")
        return False
    
    # Calculate MD5 hashes
    with open(original_file, 'rb') as f:
        original_hash = hashlib.md5(f.read()).hexdigest()
    
    with open(reconstructed_file, 'rb') as f:
        reconstructed_hash = hashlib.md5(f.read()).hexdigest()
    
    if original_hash == reconstructed_hash:
        print(f"✓ File integrity verified: {original_hash}")
        return True
    else:
        print(f"✗ File integrity mismatch:")
        print(f"  Original: {original_hash}")
        print(f"  Reconstructed: {reconstructed_hash}")
        return False

def test_basic_functionality():
    """Test basic encoding/decoding functionality."""
    print("\n=== Testing Basic Functionality ===")
    
    # Create output directory
    test_dir = os.path.join(BASE_PATH, "basic_test")
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)
    os.makedirs(test_dir, exist_ok=True)
    
    # Create a test file (5KB)
    test_file = os.path.join(test_dir, "input.bin")
    create_test_file(test_file, 5 * 1024)
    
    # Initialize MolFS interface
    molfs_dev = MolFSDev()
    molfs_dev.set_block_size(5 * 1024)  # 5KB blocks
    
    # Register custom primers - must be exactly 22 bases
    molfs_dev.register_primers(1, 0, "CTACACGACGCTCTTCCGATCT", "AGATCGGAAGAGCGGTTCAGCA")
    
    # Set the pool and block
    molfs_dev.Pool = 1
    molfs_dev.Block = 0
    
    # Encode the file
    encoded_file = os.path.join(test_dir, "encoded.csv")
    success, _ = molfs_dev.encode(test_file, encoded_file)
    
    if not success:
        print("❌ Encoding failed")
        return False
    
    # Decode the file
    decoded_file = os.path.join(test_dir, "decoded.bin")
    success, block_info = molfs_dev.decode(encoded_file, decoded_file)
    
    if not success:
        print("❌ Decoding failed")
        return False
    
    # Verify file integrity
    return verify_file_integrity(test_file, decoded_file)

def test_multi_block():
    """Test multi-block encoding/decoding."""
    print("\n=== Testing Multi-Block Functionality ===")
    
    # Create output directory
    test_dir = os.path.join(BASE_PATH, "multi_block_test")
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)
    os.makedirs(test_dir, exist_ok=True)
    
    # Create a test file (15KB - should be 3 blocks at 5KB each)
    test_file = os.path.join(test_dir, "input.bin")
    create_test_file(test_file, 15 * 1024)
    
    # Initialize MolFS interface
    molfs_dev = MolFSDev()
    molfs_dev.set_block_size(5 * 1024)  # 5KB blocks
    
    # Register unique primers for different blocks - must be exactly 22 bases
    molfs_dev.register_primers(1, 0, "CTACACGACGCTCTTCCGATCT", "AGATCGGAAGAGCGGTTCAGCA")
    molfs_dev.register_primers(1, 1, "CTACACGACACTTTTCCGATCT", "AGATCGGAAGAGCGGAACAGCA")
    molfs_dev.register_primers(1, 2, "CTACACGACGCTAATCCGATCT", "AGATCGGAAGAGCGGGTCAGCA")
    
    # Define a simple distribution strategy - all blocks to pool 1
    def simple_strategy(block_idx, total_blocks):
        return [1]  # All blocks go to Pool 1
    
    # Encode the file across blocks
    encoded_dir = os.path.join(test_dir, "encoded")
    block_distribution = molfs_dev.encode_file(test_file, encoded_dir, simple_strategy)
    
    # Get all block files
    block_files = [info["file"] for _, info in block_distribution.items()]
    
    # Reconstruct the file
    reconstructed_file = os.path.join(test_dir, "reconstructed.bin")
    success, file_info = molfs_dev.reconstruct_file(block_files, reconstructed_file)
    
    if not success:
        print("❌ Reconstruction failed")
        return False
    
    # Verify file integrity
    return verify_file_integrity(test_file, reconstructed_file)

def test_redundant_storage():
    """Test redundant block storage across different pools."""
    print("\n=== Testing Redundant Storage ===")
    
    # Create output directory
    test_dir = os.path.join(BASE_PATH, "redundant_test")
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)
    os.makedirs(test_dir, exist_ok=True)
    
    # Create a test file (15KB - should be 3 blocks at 5KB each)
    test_file = os.path.join(test_dir, "input.bin")
    create_test_file(test_file, 15 * 1024)
    
    # Initialize MolFS interface
    molfs_dev = MolFSDev()
    molfs_dev.set_block_size(5 * 1024)  # 5KB blocks
    
    # Register unique primers for different Pool/Block combinations - must be exactly 22 bases
    # Block 1 will be in both Pool 1 and Pool 2
    molfs_dev.register_primers(1, 0, "CTACACGACGCTCTTCCGATCT", "AGATCGGAAGAGCGGTTCAGCA")
    molfs_dev.register_primers(1, 1, "CTACACGACACTTTTCCGATCT", "AGATCGGAAGAGCGGAACAGCA")
    molfs_dev.register_primers(2, 1, "CTACACGACGCTTAACCGATCT", "AGATCGGAAGAGCGGATCAGCA")
    molfs_dev.register_primers(1, 2, "CTACACGACGCTAATCCGATCT", "AGATCGGAAGAGCGGGTCAGCA")
    
    # Define a distribution strategy with redundancy for block 1
    def redundant_strategy(block_idx, total_blocks):
        if block_idx == 1:
            # Block 1 goes to both Pool 1 and Pool 2
            return [1, 2]
        else:
            # Other blocks just go to Pool 1
            return [1]
    
    # Encode the file across pools
    encoded_dir = os.path.join(test_dir, "encoded")
    block_distribution = molfs_dev.encode_file(test_file, encoded_dir, redundant_strategy)
    
    # Simulate failure of Pool 1, Block 1 by excluding it
    # Get all files except Pool 1, Block 1
    block_files = [info["file"] for (pool, block), info in block_distribution.items() 
                  if not (pool == 1 and block == 1)]
    
    # Reconstruct with "damaged" Pool 1
    reconstructed_file = os.path.join(test_dir, "reconstructed.bin")
    success, file_info = molfs_dev.reconstruct_file(block_files, reconstructed_file)
    
    if not success:
        print("❌ Reconstruction failed")
        return False
    
    # Verify file integrity
    return verify_file_integrity(test_file, reconstructed_file)

def run_all_tests():
    """Run all tests and report results."""
    # Create base directory if it doesn't exist
    os.makedirs(BASE_PATH, exist_ok=True)
    
    tests = [
        ("Basic Functionality", test_basic_functionality),
        ("Multi-Block", test_multi_block),
        ("Redundant Storage", test_redundant_storage)
    ]
    
    results = {}
    
    print("=== Church Algorithm Implementation Test Suite ===")
    print(f"Using base path: {BASE_PATH}")
    print(f"Running {len(tests)} tests...\n")
    
    for name, test_func in tests:
        print(f"\n{'=' * 50}")
        print(f"Running test: {name}")
        print(f"{'=' * 50}")
        
        try:
            result = test_func()
            results[name] = result
        except Exception as e:
            print(f"❌ Test '{name}' raised an exception: {e}")
            results[name] = False
    
    # Print summary
    print("\n\n=== Test Results Summary ===")
    all_passed = True
    for name, result in results.items():
        status = "✅ PASSED" if result else "❌ FAILED"
        print(f"{status} - {name}")
        if not result:
            all_passed = False
    
    return all_passed

if __name__ == "__main__":
    success = run_all_tests()
    
    if success:
        print("\n🎉 All tests passed! The implementation is working correctly.")
        sys.exit(0)
    else:
        print("\n❌ Some tests failed. Please check the logs above for details.")
        sys.exit(1)
