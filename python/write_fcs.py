#!/usr/bin/python
# -*- coding: utf-8 -*-
# This file was modified from an existing repository and 99% of props go the the author of
# https://github.com/ZELLMECHANIK-DRESDEN/fcswrite/blob/master/fcswrite/fcswrite.py
"""Write .fcs files for flow cytometry"""

import numpy as np
import pathlib
import struct
import warnings

def write_fcs(filename, chn_names, data,
              chn_names_long=None,
              text_kw_pr={},
              endianness=None,
              compat_chn_names=False,
              compat_copy=False,
              compat_negative=False,
              compat_percent=False,
              compat_max_int16=10000):
    """Write numpy data to an .fcs file (FCS3.0 file format)
    Parameters
    ----------
    filename: str or pathlib.Path
        Path to the output .fcs file
    chn_names: list of str, length C
        Names of the output channels
    data: 2d ndarray of shape (N,C)
        The numpy array data to store as .fcs file format.
    chn_names_long: list of str, length C
        Long Names of the output channels
    text_kw_pr: dict
        User-defined, optional key-value pairs that are stored
        in the primary TEXT segment
    endianness: str|None
        Set to "little" or "big" to explicitly define the byte 
        order used. If None, the endianness is inherited from the
        $BYTEORD key in text_kw_pr
    compat_chn_names: bool
        Compatibility mode for 3rd party flow analysis software:
        The characters " ", "?", and "_" are removed in the output
        channel names.
    compat_copy: bool
        Do not override the input array `data` when modified in
        compatibility mode.
    compat_negative: bool
        Compatibliity mode for 3rd party flow analysis software:
        Flip the sign of `data` if its mean is smaller than zero.
    compat_percent: bool
        Compatibliity mode for 3rd party flow analysis software:
        If a column in `data` contains values only between 0 and 1,
        they are multiplied by 100.
    compat_max_int16: int
        Compatibliity mode for 3rd party flow analysis software:
        If a column in `data` has a maximum above this value,
        then the display-maximum is set to 2**15.
    Notes
    -----
    - These commonly used unicode characters are replaced: "µ", "²"
    - If the input data contain NaN values, the corresponding rows
      are excluded due to incompatibility with the FCS file format.
    """
    # Put this in a fresh dict since we modify it later in the function
    _text_kw_pr = text_kw_pr.copy()
    # Drop the keys that will need to be filled at write-time by this function
    for k in ['__header__', 
              '$BEGINANALYSIS', '$ENDANALYSIS', '$BEGINSTEXT', '$ENDSTEXT', 
              '$BEGINDATA', '$ENDDATA', '$DATATYPE', '$MODE', 
              '$NEXTDATA', '$TOT', '$PAR']:
        _text_kw_pr.pop(k, None)

    filename = pathlib.Path(filename)
    if not isinstance(data, np.ndarray):
        data = np.array(data, dtype=float)
    # remove rows with nan values
    nanrows = np.isnan(data).any(axis=1)
    if np.sum(nanrows):
        msg = "Rows containing NaNs are not written to {}!".format(filename)
        warnings.warn(msg)
        data = data[~nanrows]

    if endianness == "little":
        # use little endian
        byteord = '1,2,3,4'
    elif endianness == "big":
        # use big endian
        byteord = '4,3,2,1'
    else:
        try:
            byteord = _text_kw_pr.pop('$BYTEORD')
        except KeyError:
            raise ValueError('Cannot have `endianness=None` if `$BYTEORD` is not a key in '
                               '`text_kw_pr`')
        else:
            if byteord not in ('4,3,2,1', '1,2,3,4'):
                raise ValueError('`text_kw_pr["$BYTEORD"]` can only be one of `1,2,3,4` or '
                                 '`4,3,2,1`.')

    msg = "length of `chn_names` must match length of 2nd axis of `data`"
    assert len(chn_names) == data.shape[1], msg

    rpl = [["µ", "u"],
           ["²", "2"],
           ]

    if compat_chn_names:
        # Compatibility mode: Clean up headers.
        rpl += [[" ", ""],
                ["?", ""],
                ["_", ""],
                ]

    for ii in range(len(chn_names)):
        for (a, b) in rpl:
            chn_names[ii] = chn_names[ii].replace(a, b)
            if chn_names_long is not None:
                chn_names_long[ii] = chn_names_long[ii].replace(a, b)

    # Data with values between 0 and 1
    pcnt_cands = np.apply_along_axis(lambda x: (x.min() >= 0 and x.max() <= 1), 0, data)
    pcnt_cands = [i for i, p in enumerate(pcnt_cands) if p]
    if compat_percent and pcnt_cands:
        # Compatibility mode: Scale values b/w 0 and 1 to percent
        if compat_copy:
            # copy if requested
            data = data.copy()
        for ch in pcnt_cands:
            data[:, ch] *= 100

    if compat_negative:
        toflip = np.apply_along_axis(lambda x: np.mean(x)  < 0, 0, data)
        toflip = [i for i, v in enumerate(toflip) if v]
        if len(toflip):
            if compat_copy:
                # copy if requested
                data = data.copy()
            for ch in toflip:
                data[:, ch] *= -1

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
        # Set data maximum to that of int16
        if (compat_max_int16 and
            np.max(data[:, jj]) > compat_max_int16 and
                np.max(data[:, jj]) < 2**15):
            pnrange = int(2**15)
        # Set range for data with values between 0 and 1
        elif jj in pcnt_cands:
            if compat_percent:  # scaled to 100%
                pnrange = 100
            else:  # not scaled
                pnrange = 1
        # default: set range to maxium value found in column
        else:
            pnrange = int(abs(np.max(data[:, jj])))
        # TODO:
        # - Set log/lin
        # Using pop on _text_kw_pr will remove it from the dict if it exists else it will use a 
        # default value . 
        data_str = '/'.join(['',
                            '$P{0}B'.format(jj+1),
                            str(_text_kw_pr.pop('$P{0}B'.format(jj+1), '32')),
                            '$P{0}E'.format(jj+1),
                            str(_text_kw_pr.pop('$P{0}E'.format(jj+1), '0,0')),
                            '$P{0}N'.format(jj+1),
                            str(_text_kw_pr.pop('$P{0}N'.format(jj+1), chn_names[jj]))])
        if chn_names_long is not None:
            data_str = '/'.join([data_str,
                                '$P{0}S'.format(jj+1),
                                str(_text_kw_pr.pop('$P{0}S'.format(jj+1), chn_names_long[jj]))])
        data_str = '/'.join([data_str,
                            '$P{0}R'.format(jj+1),
                            str(_text_kw_pr.pop('$P{0}R'.format(jj+1), pnrange)),
                            '$P{0}D'.format(jj+1),
                            str(_text_kw_pr.pop('$P{0}D'.format(jj+1), 'Linear')),
                            '$P{0}G'.format(jj+1),
                            str(_text_kw_pr.pop('$P{0}G'.format(jj+1), '1')),
                            ])
        TEXT += data_str

    # Finally, add any remaining, additional key-value pairs provided by the user.
    for key in sorted(_text_kw_pr.keys()):
        TEXT += '/{0}/{1}'.format(key, text_kw_pr[key])

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
    with filename.open("wb") as fd:
        fd.write(HEADER.encode("ascii", "replace"))
        fd.write(TEXT.encode("ascii", "replace"))
        fd.write(DATA)
        fd.write(b'00000000')