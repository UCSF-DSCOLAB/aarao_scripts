#!/bin/env python3
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--prefixB', action='store_true', help='Prefix -B to each folder? Useful for singularity')
    parser.add_argument('folders', nargs='+', help='list of folders to process.')
    params = parser.parse_args()

    folders = params.folders
    # Sort longest path to shortest
    folders = sorted(folders, key=lambda x: -x.count('/'))
    out_folders = []
    
    prefix = '-B ' if params.prefixB else ''
    for i in range(len(folders)):
        for j in range(i+1, len(folders)):
            if folders[i].startswith(folders[j]):
                break
        else:
            out_folders.append(f'{prefix}{folders[i]}')

    print(*out_folders, sep=' ')
if __name__ == '__main__':
    main()
