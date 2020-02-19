import argparse
import gzip
import os
import pysam


def generate_fasta(bamfile, outfile, noCB_outfile, update_verbosity=10000):
    udict = {}
    sam_handle = pysam.Samfile(bamfile, 'rb')

    fa_handle = gzip.open(outfile, 'wt')
    noCBfa_handle = gzip.open(noCB_outfile, 'wt')

    noCB_readnum = bad_umi = umi_dup = good_read = num_cells = 0
    try:
        i=0
        for read in sam_handle.fetch(until_eof=True):
            i += 1
            if i % update_verbosity == 0:
                print('Processed {} reads'.format(i))

            if not read.is_unmapped:
                continue
            rdict = dict(read.tags)
            if 'fr' in rdict:
                # This is a feature barcoding read and thus is useless
                print('Found a Feature barcode')
                continue
            elif 'CB' not in rdict:
                out_read = ''.join([(b if ord(q)>=25 else 'N')
                                    for b, q in zip(read.seq, read.qual)])
                noCBfa_handle.write('@no_cell_barcode_read{}\n{}\n'.format(noCB_readnum, out_read))
                noCB_readnum += 1
            else:
                if 'UB' not in rdict:
                    # Bad Umi... skip
                    bad_umi += 1
                    continue
                if rdict['CB'] in udict:
                    if rdict['UB'] in udict[rdict['CB']]:
                        # This is a umi duplicate
                        print('Found a UMI duplicate')
                        umi_dup += 1
                        continue
                    else:
                        good_read += 1
                        udict[rdict['CB']].add(rdict['UB'])
                else:
                    good_read += 1
                    num_cells += 1
                    udict[rdict['CB']] = {rdict['UB']}

                out_read = ''.join([(b if ord(q)>=25 else 'N')
                                    for b, q in zip(read.seq, read.qual)])
                fa_handle.write('@{}:{}\n{}\n'.format(rdict['CB'],
                                                    rdict['UB'],
                                                    out_read))
    except:
        print(read.to_dict())
        print("######")
        print(read.tags)
        raise
    finally:
        sam_handle.close()
        fa_handle.close()
        noCBfa_handle.close()
        print('Identified:',
              'Good Cells:{}'.format(num_cells),
              'Good reads:{}'.format(good_read),
              'Cells without proper barcodes: {}'.format(noCB_readnum),
              'Reads with bad UMIs:{}'.format(bad_umi),
              'UMI duplicates:{}'.format(umi_dup),
              sep='\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('bamfile')
    parser.add_argument('outfile_prefix')
    parser.add_argument('--update_verbosity', type=int, default=10000)
    params = parser.parse_args()

    assert os.path.exists(params.bamfile)
    assert not os.path.isdir(params.outfile_prefix)
    outfile = ''.join([params.outfile_prefix, '.fasta.gz'])
    assert not os.path.exists(outfile)
    noCB_outfile = ''.join([params.outfile_prefix, '_noCB.fasta.gz'])
    assert not os.path.exists(noCB_outfile)

    generate_fasta(params.bamfile, outfile, noCB_outfile, params.update_verbosity)

if __name__ == '__main__':
    main()