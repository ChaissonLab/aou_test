#!/usr/bin/env python
import sys
tsv = open(sys.argv[1])
n=int(sys.argv[2])
start=0
if (len(sys.argv) > 3):
    start = int(sys.argv[3])
entries=[]
q="\""
i=0
tsv.readline()
for line in tsv:
    vals=line.strip().split(",")
    if i >= start:
        entries.append(vals)
    i+=1
    if i == n+start:
        break

print("""
{
  "RunCtyperBatch.input_bams": [ """)
for i in range(0,len(entries)-1):
    print(q + entries[i][1] + q+",")
print(q + entries[-1][1] + q + "],")

print(""""RunCtyperBatch.input_bai_files": [""")
for i in range(0,len(entries)-1):
    print(q + entries[i][2] + q+",")
print(q + entries[-1][2] + q + "],")

print("""    "RunCtyperBatch.output_tab_files": [""")
for i in range(0,len(entries)-1):
    print(q + entries[i][0] + ".tab" + q +",")
print(q + entries[-1][0] + ".tab" + q + "],")

print(""""RunCtyperBatch.matrix": "gs://fc-secure-eb06d0c9-fe8f-427b-8acf-caac7e60b975/ALL_kmatrix01.txt",""")
print(""""RunCtyperBatch.matrixIndex": "gs://fc-secure-eb06d0c9-fe8f-427b-8acf-caac7e60b975/ALL_kmatrix01.txt.index",""")
print(""""RunCtyperBatch.matrixBin": "gs://fc-secure-eb06d0c9-fe8f-427b-8acf-caac7e60b975/ALL_kmatrix01.txt.bin"\n}""")

