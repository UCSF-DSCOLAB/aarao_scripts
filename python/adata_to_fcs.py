#!/usr/bin/python
# -*- coding: utf-8 -*-
# This file was modified from an existing repository and 99% of props go the the author of
# https://github.com/ZELLMECHANIK-DRESDEN/fcswrite/blob/master/fcswrite/fcswrite.py
"""Write .fcs files for flow cytometry"""

import numpy as np
import os
import pandas as pd
import struct
import warnings

from datetime import datetime
from scipy.sparse.csr import csr_matrix


def fcs_from_adata(filename, 
                   adata,
                   raw=False,
                   scaled_col=None,
                   embeddings=None,
                   chn_names_short_col=None,
                   event_metadata_cols=None,
                   event_metadata_prefix='',
                   add_cellular_barcode=False,
                   cellular_barcode_additive=0,
                   cellular_barcode_names=None,
                   sample_name_col=None,
                   text_kw_pr=None):
    """Write numpy data to an .fcs file (FCS3.0 file format)
    Parameters
    ----------
    filename: str or pathlib.Path
        Path to the output .fcs file
    adata: 2d ndarray of shape (N,C)
        The numpy array data to store as .fcs file format.
    raw: bool
        Use the raw data from the adata object?
    scaled_col: 
        If not raw, a boolean column in adata.var that describes columns that should be Log scaled (True)
    embeddings: str or list of str
        List of dimensionality reductions to append into the fcs. Will only use the first 2 dims.
    chn_names_short_col: str
        Column name in adata.var that has the short name for the channel
    event_metadata_cols: str or list(str)
        Column name(s) in adata.obs that have event-level metadata.
    event_metadata_prefix: str or list(str)
        Prefix to append to the event labels. Must be the same length/type as event_metadata_cols
    add_cellular_barcode: bool
        Should we add a channel with a unique numeric cellular barcode (Useful if pulling the data 
        back into python later from cellengine, say)
    cellular_barcode_additive: int
        If we are adding cellular barcodes, should we add a value to the barcode? This can be useful
        to inject library information into the barcode (prior to merging, say)
    cellular_barcode_names: tuple(str, str)
        Short and long names respectively for the cellular barcode channel
    sample_name_col: str
        Column name in adata.obs that has the sample annotaitons per event
    text_kw_pr: dict
        User-defined, optional key-value pairs that are stored
        in the primary TEXT segment
    Notes
    -----
    - These commonly used unicode characters are replaced: "µ", "²"
    - If the input data contain NaN values, the corresponding rows
      are excluded due to incompatibility with the FCS file format.
    """
    if add_cellular_barcode:
        assert isinstance(cellular_barcode_additive, int)
        if cellular_barcode_names is not None:
            assert len(cellular_barcode_names) == 2
            assert all([isinstance(cbn, str) for cbn in cellular_barcode_names])
        else:
            cellular_barcode_names = ['barcode', 'barcode']
    chn_names = adata.var.index.tolist()
    chn_names_original = adata.var.index.tolist()
    if chn_names_short_col is None:
        chn_names_short = chn_names
    else:
        chn_names_short = adata.var[chn_names_short_col].tolist()
    if text_kw_pr is None:
        text_kw_pr = {}
    
    if raw or scaled_col is None:
        scale = ['Linear'] * len(chn_names)
    elif isinstance(scaled_col, str):
        # This may not work :shrug:
        assert scaled_col in adata.var.columns and pd.api.types.is_bool_dtype(adata.var[scaled_col])
        scale = ['Logarithmic' if v else "Linear" for v in adata.var[scaled_col]]
    if embeddings is None:
        embeddings = []
    elif isinstance(embeddings, str):
        embeddings = [embeddings]
    elif isinstance(embeddings, list):
        pass
    else:
        raise RuntimeError('embeddings must be of type `None`, `str`, or `list(str)`')

    assert all([f'X_{dr}' in adata.obsm for dr in embeddings])

    # Will create event_metadata as a dict of column: prefix_for_column
    if event_metadata_cols is None:
        event_metadata = {}
    elif isinstance(event_metadata_cols, str):
        assert isinstance(event_metadata_prefix, str), 'event_metadata_prefix and event_metadata_cols must be of the same type'
        event_metadata = {event_metadata_cols: event_metadata_prefix}
    elif isinstance(event_metadata_cols, list):
        assert len(event_metadata_cols) == len(event_metadata_prefix), 'event_metadata_prefix and event_metadata_cols must be of the same length'
        assert all([isinstance(x, str) for x in event_metadata_cols]), 'event_metadata_cols must be all strings'
        assert all([isinstance(x, str) for x in event_metadata_prefix]), 'event_metadata_prefix must be all strings'
        event_metadata = dict(zip(event_metadata_cols, event_metadata_prefix))
    else:
        raise RuntimeError('event_metadata_cols must be of type `None`, `str`, or `list(str)`')
    
    if sample_name_col not in event_metadata:
        event_metadata[sample_name_col] = ''
    
    # Put this in a fresh dict since we modify it later in the function
    _text_kw_pr = text_kw_pr.copy()
    # Drop the keys that will need to be filled at write-time by this function
    for k in ['__header__', 
              '$BEGINANALYSIS', '$ENDANALYSIS', '$BEGINSTEXT', '$ENDSTEXT', 
              '$BEGINDATA', '$ENDDATA', '$DATATYPE', '$MODE', 
              '$NEXTDATA', '$TOT', '$PAR', '$FIL']:
        _text_kw_pr.pop(k, None)

    filename = os.path.abspath(filename)
    if raw:
        if isinstance(adata.raw.X, csr_matrix):
            data = np.array(adata.raw.X.todense(), dtype=float)
        else:
            data = np.array(adata.raw.X, dtype=float)
    else:
        if isinstance(adata.X, csr_matrix):
            data = np.array(adata.X.todense(), dtype=float)
        else:
            data = np.array(adata.X, dtype=float)
    #TODO: These 3 can just be one large array to reduce redundancy
    # Add the metadata to the data at this point
    metadata_arr = np.empty((data.shape[0],0), np.float64)
    metadata_labels = []
    metadata_key = {}
    for label, prefix in event_metadata.items():
        metadata_labels.append(label)
        if pd.api.types.is_integer_dtype(adata.obs[label]):
             metadata_arr = np.column_stack((metadata_arr, adata.obs[label].to_numpy()))
             metadata_key[label] = {x: f'{prefix}{x}' for x in adata.obs[label].unique()}
        elif pd.api.types.is_categorical_dtype(adata.obs[label]):
            metadata_arr = np.column_stack((metadata_arr, adata.obs[label].cat.codes))
            metadata_key[label] = {i: f'{prefix}{x}' for i, x in enumerate(adata.obs[label].cat.categories)} 
        else:
            # cast to categiorical and work from there.
            temp = adata.obs[label].astype('category').copy()
            metadata_arr = np.column_stack((metadata_arr, temp.cat.codes))
            metadata_key[label] = {i: f'{prefix}{x}' for i, x in enumerate(temp.cat.categories)}                
    if metadata_arr.shape[1] > 0:
        data = np.hstack((data, metadata_arr))
        chn_names.extend(metadata_labels)
        chn_names_short.extend(metadata_labels)
        scale.extend(['Linear'] * len(metadata_labels))
    # Make some space since this may be big
    del(metadata_arr)
    # Now add the cellular barcode, if any.
    if add_cellular_barcode:
        cb_arr = np.array(cellular_barcode_additive + np.arange(start=1, stop=data.shape[0]+1, 
                                                                step=1, dtype=np.int64))
        data = np.hstack((data, cb_arr.reshape(cb_arr.shape[0], 1)))
        chn_names.append(cellular_barcode_names[1])
        chn_names_short.append(cellular_barcode_names[0])
        scale.append('Linear')
        # Delete unnecessary array
        del(cb_arr)
    # Add the dimensionality reductions if any
    dimred_arr = np.empty((data.shape[0],0), np.float64)
    dimred_labels = []
    for dr in embeddings:
        dimred_labels.extend([f'{dr}_1', f'{dr}_2'])
        dimred_arr = np.hstack((dimred_arr, adata.obsm[f'X_{dr}'][:,:2]))
    if dimred_arr.shape[1] > 0:
        data = np.hstack((data, dimred_arr))
        chn_names.extend(dimred_labels)
        chn_names_short.extend(dimred_labels)
        scale.extend(['Linear'] * len(dimred_labels))
    # Make some space since this may be big
    del(dimred_arr)

    # remove rows with nan values
    nanrows = np.isnan(data).any(axis=1)
    if np.sum(nanrows):
        msg = "Rows containing NaNs are not written to {}!".format(filename)
        warnings.warn(msg)
        data = data[~nanrows]

    # Always use big endianness
    _text_kw_pr.pop('$BYTEORD', None)
    byteord = '4,3,2,1'

    msg = f"length of `chn_names` ({len(chn_names)}) must match length of 2nd axis of `data`({data.shape[1]})"
    assert len(chn_names) == data.shape[1], msg

    rpl = [["µ", "u"],
           ["²", "2"],
           ]

    for ii in range(len(chn_names)):
        for (a, b) in rpl:
            chn_names[ii] = chn_names[ii].replace(a, b)
            chn_names_short[ii] = chn_names_short[ii].replace(a, b)

    # Data with values between 0 and 1
    pcnt_cands = np.apply_along_axis(lambda x: (x.min() >= 0 and x.max() <= 1), 0, data)
    pcnt_cands = [i for i, p in enumerate(pcnt_cands) if p and chn_names[i] in chn_names_original]

    # DATA segment
    data1 = data.flatten().tolist()
    DATA = struct.pack(f'>{data.shape[0] * data.shape[1]}f', *data.flatten())

    # TEXT segment
    header_size = 256

    TEXT = '/$BEGINANALYSIS/0/$ENDANALYSIS/0'
    TEXT += '/$BEGINSTEXT/0/$ENDSTEXT/0'
    # Add placeholders for $BEGINDATA and $ENDDATA, because we don't
    # know yet how long TEXT is.
    TEXT += '/$BEGINDATA/{data_start_byte}/$ENDDATA/{data_end_byte}'
    TEXT += '/$BYTEORD/{0}/$DATATYPE/F'.format(byteord)
    TEXT += '/$MODE/L/$NEXTDATA/0/$TOT/{0}'.format(data.shape[0])
    TEXT += '/$PAR/{0}'.format(data.shape[1])

    # Check for content of data columns and set range
    for jj in range(data.shape[1]):
        # Set range for data with values between 0 and 1
        if jj in pcnt_cands:
            pnrange = 1
        # default: set range to maxium value found in column
        else:
            pnrange = int(abs(np.max(data[:, jj])))
        # Using pop on _text_kw_pr will remove it from the dict if it exists else it will use a 
        # default value . 
        data_str = '/'.join(['',
                            '$P{0}B'.format(jj+1),
                            str(_text_kw_pr.pop('$P{0}B'.format(jj+1), '32')),
                            '$P{0}E'.format(jj+1),
                            str(_text_kw_pr.pop('$P{0}E'.format(jj+1), '0,0')),
                            '$P{0}N'.format(jj+1),
                            str(_text_kw_pr.pop('$P{0}N'.format(jj+1), chn_names[jj])),
                            '$P{0}S'.format(jj+1),
                            str(_text_kw_pr.pop('$P{0}S'.format(jj+1), chn_names_short[jj])),
                            '$P{0}R'.format(jj+1),
                            str(_text_kw_pr.pop('$P{0}R'.format(jj+1), pnrange)),
                            '$P{0}D'.format(jj+1),
                            str(_text_kw_pr.pop('$P{0}D'.format(jj+1), scale[jj])),
                            '$P{0}G'.format(jj+1),
                            str(_text_kw_pr.pop('$P{0}G'.format(jj+1), '1')),
                            ])
        TEXT += data_str

    # Finally, add any remaining, additional key-value pairs provided by the user.
    # First add some fixed values to this dict
    _text_kw_pr['$FIL'] = os.path.basename(filename)
    if '$EXP' not in _text_kw_pr:
        _text_kw_pr['$EXP'] = os.environ["USER"]
    _text_kw_pr['$COM'] = [
        f'This file was generated on {datetime.now().strftime("%Y-%m-%d at %H:%M:%S")} using an automated script.',
        f'Custom Dimensionality Reduction [i.e. Embeddding] parameters included: {embeddings}.',
        f'Custom metadata containing parameters included: {metadata_labels}.',
        'The mapping of metadata parameters is as follows:']

    for k, v in metadata_key.items():
        _text_kw_pr['$COM'].append(f'{k} :: {str(v).replace("{", "[").replace("}", "]")}. ')

    _text_kw_pr['$COM'] = '\n '.join(_text_kw_pr['$COM'])
    
    for key in sorted(_text_kw_pr.keys()):
        TEXT += '/{0}/{1}'.format(key, _text_kw_pr[key])

    #Add a trailing slash
    TEXT += '/'

    # SET $BEGINDATA and $ENDDATA using the current size of TEXT plus padding.
    text_padding = 47  # for visual separation and safety
    data_start_byte = header_size + len(TEXT) + text_padding
    data_end_byte = data_start_byte + len(DATA) - 1
    TEXT = TEXT.format(data_start_byte=data_start_byte,
                       data_end_byte=data_end_byte)
    lentxt = len(TEXT)
    # Pad TEXT segment with spaces until data_start_byte
    TEXT = TEXT.ljust(data_start_byte - header_size, " ")
    # HEADER segment
    ver = 'FCS3.0'
    # Get values for the analysis segment
    textfirst = '{0: >8}'.format(header_size)
    textlast = '{0: >8}'.format(lentxt + header_size - 1)
    # Starting with FCS 3.0, data segment can end beyond byte 99,999,999,
    # in which case a zero is written in each of the two header fields (the
    # values are given in the text segment keywords $BEGINDATA and $ENDDATA)
    if data_end_byte <= 99999999:
        datafirst = '{0: >8}'.format(data_start_byte)
        datalast = '{0: >8}'.format(data_end_byte)
    else:
        datafirst = '{0: >8}'.format(0)
        datalast = '{0: >8}'.format(0)
    # Spoof values for the analysis segment
    anafirst = '{0: >8}'.format(0)
    analast = '{0: >8}'.format(0)
    # Generate the header segment
    HEADER = '{0: <256}'.format(ver + '    '
                                + textfirst
                                + textlast
                                + datafirst
                                + datalast
                                + anafirst
                                + analast)
    # Write data
    with open(filename, "wb") as fd:
        fd.write(HEADER.encode("ascii", "replace"))
        fd.write(TEXT.encode("ascii", "replace"))
        fd.write(DATA)
        fd.write(b'00000000')