import argparse
import gzip
import os
import pysam


def generate_fasta(bamfile, outfile, update_verbosity=10000):
    udict = {}
    sam_handle = pysam.Samfile(bamfile, 'rb')

    noCB_readnum = 0
    bad_umi = 0
    umi_dup = 0
    good_read = 0
    num_cells = 0
    unmapped = 0
    duplicate = 0
    feature_barcode = 0
    low_map_qual = 0

    reads_processed = 0
    try:
        for start, stop in [(1, 10001), (2781479,56887903)]:
            #i = 0
            for read in sam_handle.fetch('Y', start=start, stop=stop):
                reads_processed += 1
                #if i == 5000:
                #    break
                if reads_processed % update_verbosity == 0:
                    print('Processed {} reads'.format(reads_processed))

                if read.is_unmapped:
                    unmapped += 1
                    continue
                if read.is_duplicate:
                    duplicate += 1
                    continue
                if read.mapping_quality < 255:
                    low_map_qual += 1
                    continue
                #i += 1

                rdict = dict(read.tags)
                if 'fr' in rdict:
                    feature_barcode += 1
                    # This is a feature barcoding read and thus is useless
                    continue
                elif 'CB' not in rdict:
                    noCB_readnum += 1
                else:
                    if 'UB' not in rdict:
                        # Bad Umi... skip
                        bad_umi += 1
                        continue
                    if rdict['CB'] in udict:
                        if rdict['UB'] in udict[rdict['CB']]:
                            # This is a umi duplicate
                            umi_dup += 1
                            continue
                        else:
                            #print(read.mapping_quality)
                            #print(read.is_secondary)
                            #print(read.is_supplementary)
                            good_read += 1
                            udict[rdict['CB']].add(rdict['UB'])
                    else:
                        #print(read.mapping_quality)
                        #print(read.is_secondary)
                        #print(read.is_supplementary)
                        good_read += 1
                        num_cells += 1
                        udict[rdict['CB']] = {rdict['UB']}
    except:
        print('Failed on read {}'.format(read.query_name))
        print(read.to_dict())
        print('######')
        print(read.tags)
        raise
    finally:
        sam_handle.close()
        print('Identified:',              
              'Good Cells:{}'.format(num_cells),
              'Cells without proper barcodes: {}'.format(noCB_readnum),
              'Total reads:{}'.format(reads_processed),
              'Good reads:{}'.format(good_read),
              'Reads with bad UMIs:{}'.format(bad_umi),
              'UMI duplicates:{}'.format(umi_dup), 
              'Unmapped reads:{}'.format(unmapped),
              'Duplicated reads:{}'.format(duplicate),
              'Reads identified as feature barcodes:{}'.format(feature_barcode),
              'Reads with low map qual:{}'.format(low_map_qual),
              sep='\n')
    cdict = {x: len(y) for x,y in udict.items()}
    
    with open(outfile, 'w') as off:
        print('barcode', 'ChrY.NonPARY.count', sep='\t', file=off)
        for bc, cnt in cdict.items():
            print(bc, cnt, sep='\t', file=off)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bamfile')
    parser.add_argument('outfile')
    parser.add_argument('--overwrite', action="store_true")
    parser.add_argument('--update_verbosity', type=int, default=10000)
    params = parser.parse_args()

    bamfile = os.path.abspath(params.bamfile)
    assert os.path.exists(bamfile)
    
    outfile = os.path.abspath(params.outfile)
    if os.path.exists(params.outfile):
        if params.overwrite:
            print('Overwriting existing file: {}'.format(outfile))

        else:
            raise RuntimeError('File exists ({}). specify --overwrite to overwrite '.format(outfile))

    generate_fasta(bamfile, outfile, params.update_verbosity)

if __name__ == '__main__':
    main()