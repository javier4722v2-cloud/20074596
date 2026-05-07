import gzip
import os
import csv
import sys
import math
from collections import defaultdict
from scipy import stats as scipy_stats

base = r'C:\Users\javie\.openclaw\workspace\lumus\data\gse81622'
gzfile = os.path.join(base, 'GSE81622_non-normalized.txt.gz')
txtfile = os.path.join(base, 'GSE81622_non-normalized.txt')

# Step 1: Decompress
print("Decompressing...")
with gzip.open(gzfile, 'rb') as f_in:
    data = f_in.read()
with open(txtfile, 'wb') as f_out:
    f_out.write(data)
print(f"Decompressed: {len(data)} bytes")

# Step 2: Read header to understand column structure
print("\nReading header...")
with open(txtfile, 'r', encoding='latin-1') as f:
    header = f.readline().strip()
    cols = header.split('\t')
    
# Print first few columns to understand structure
print(f"Total columns: {len(cols)}")
print(f"First 10 columns: {cols[:10]}")
print(f"Column at index 0: '{cols[0]}'")
print(f"Column at index 1: '{cols[1]}'")
print(f"Column at index 2: '{cols[2]}'")

# Look at a few data rows
print("\nSample data rows (first 3):")
with open(txtfile, 'r', encoding='latin-1') as f:
    f.readline()  # skip header
    for i, line in enumerate(f):
        if i >= 3:
            break
        fields = line.strip().split('\t')
        print(f"Row {i}: {fields[:6]}")
