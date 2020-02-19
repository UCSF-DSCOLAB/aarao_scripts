"""

Description:

A script to convert a qsub template into a usable qsub script.

    python %(prog)s tags <XXX_template.sh>

        Will identify tags for a script.

    python %(prog)s generate <XXX_template.sh> TAG1=value1 TAG2=value2 ....

        Will generate the usable script (XXX.sh) using the given tags.

Tags in the script are specified in the template using 
        `<__TAG__>` for a required tag with name "TAG"
        `<_#GAT#_>` for an optional tag with name "GAT"
and will be replaced with the value provided when using this script in `generate`
mode.

Optional tag values can be provided in 3 ways
    1. Global: In the `~/.qsubt_defaults` file that contains one "TAG\\tVALUE" pair
               per line.
    2. Local : In the template script header using the format
               "#@ TAG VALUE".
    3. CLI   : As an argument to `generate` on the Command line.

The order of preference for optional values are CLI > Local > Global
"""
from __future__ import print_function

import argparse
import os
import re

req_token_regex = r"\<__(?P<tag>([A-Z]*))__\>"
opt_token_regex = r"\<_#(?P<tag>([A-Z]*))#_\>"
token_regex = r"\<_[#_](?P<tag>([A-Z]*))[#_]_\>"

def get_tags(template_file, print_results=False):
    req_tokens = set()
    opt_tokens = set()
    opt_tokenvals = {}
    if os.path.exists(os.path.expanduser('~/.qsubt_defaults')):
        with open(os.path.expanduser('~/.qsubt_defaults')) as iff:
            for x in iff:
                x = x.strip().split()
                opt_tokenvals.update({x[0]: ' '.join(x[1:])})

    with open(template_file)as iff:
        for line in iff:
            if line.startswith('#@'):
                line = line.strip().split()
                opt_tokenvals[line[1]] = ' '.join(line[2:])
                continue
            req_tokens.update({x.group('tag') for x in re.finditer(req_token_regex, line)})
            opt_tokens.update({x.group('tag') for x in re.finditer(opt_token_regex, line)})
        if req_tokens.intersection(opt_tokens):
            tokens = ','.join(req_tokens.intersection(opt_tokens))
            raise RuntimeError("Cannot have a token(s) that is both required and optional (%s)"
                               % tokens)
        if set(opt_tokenvals) != opt_tokens:
            tokens = ','.join(opt_tokens.difference(set(opt_tokenvals)))
            raise RuntimeError("The following optional tokens do not have default values: %s"
                               % tokens)
    if print_results:
        print('The required tokens for %s are:' % template_file)
        for token in req_tokens:
            print(token)
        if opt_tokens:
            padding = max(len(x) for x in opt_tokens) + 2
            if padding < 7:
                padding = 7
            outstring = '{:%s}{}' % padding
            print('\nThe optional tokens for %s are:' % template_file)
            print(outstring.format('token', 'default_value'))
            print('-'*(padding + 13))  # padding + default_value
            for token, val in opt_tokenvals.items():
                print(outstring.format(token, val))
    return req_tokens, opt_tokenvals

def generate_qsub_file(args):
    template_file = args.pop(0)
    tag_values = {}
    for tag_value in args:
        tag, value = tag_value.split('=', 1)  # Make only 1 split
        tag_values[tag] = value

    required_tags, optional_values = get_tags(template_file)
    diff = required_tags - set(tag_values)
    if diff:
        raise RuntimeError("values for %s were not provided" % str(diff))
    infile = '\n'.join(x.rstrip() for x in open(template_file) if not x.startswith('#@'))

    # Opdate optional_values so optional tags specified in the CLI are given preference
    tag_values = dict(optional_values, **tag_values)

    while True:
        tags = {x.group('tag') for x in re.finditer(token_regex, infile)}
        if not tags:
            break
        for tag in tags:
            regex = re.compile(r"\<_[_#]" + tag + "[_#]_\>")
            infile = re.sub(regex, tag_values[tag], infile)

    with open(re.sub('_template', '', template_file), 'w') as off:
        print(infile, file=off)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('action', help='Action to conduct.', type=str, choices=['generate', 'tags'])
    
    params, extras = parser.parse_known_args()


    if len(extras) == 0:
        raise RuntimeError('Cannot process without a template file.')
    if not os.path.exists(extras[0]):
        raise RuntimeError('Provided template file does not exist.')
    if not extras[0].endswith('_template.sh'):
        raise RuntimeError('Is the input file in the form `XXX_template.sh` ?')

    if params.action == 'tags':
        get_tags(extras[0], print_results=True)
    else:
        generate_qsub_file(extras)

if __name__ == '__main__':
    main()