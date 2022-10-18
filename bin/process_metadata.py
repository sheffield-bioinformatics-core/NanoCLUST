#!/usr/bin/env python

import argparse
import shutil
import pandas as pd
import os.path
import csv

def process_metadata(file):
    meta=pd.read_excel(file, usecols=range(0,6))
    barcodes=[]
    for index, row in meta.iterrows():
        if row['Status'] not in ['discontinued', 'positive control', 'negative control']:
            barcodes.append(row['Barcode'])

    if 'positive control' in set(meta['Status']):
        print("true")
        positive_ctrl='rel_abundance_'+ meta.loc[meta['Status'] == 'positive control', 'Barcode'].item() + '_S.csv'
    else:
        positive_ctrl="none"

    if 'negative control' in set(meta['Status']):
        print("true")
        negative_ctrl='rel_abundance_'+ meta.loc[meta['Status'] == 'negative control', 'Barcode'].item() + '_S.csv'
    else: 
        negative_ctrl="none"

    print(positive_ctrl)
    print(negative_ctrl)
    return positive_ctrl, negative_ctrl, barcodes

def output_meta(pos, neg):
    if os.path.exists(pos):
        shutil.copyfile(pos, "positive_control.csv")
    else:
        with open('positive_control.csv', 'w') as file:
            pass
        print("positive control not provided")

    if os.path.exists(neg):
        shutil.copyfile(neg, "negative_control.csv")
    else:
        with open('negative_control.csv', 'w') as file:
            pass
        print("negative control not provided")

    return None

def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--metatable", required=True,
        help="Table with information on positive and negative controls")
    args = parser.parse_args()

    positive,negative,barcodes=process_metadata(args.metatable)

    barcodes_rows=zip(barcodes)
    with open("barcodes.csv", "w") as f:
        writer = csv.writer(f)
        for row in barcodes_rows:
            writer.writerow(row)

    output_meta(positive, negative)
    

if __name__ == "__main__":
    main()