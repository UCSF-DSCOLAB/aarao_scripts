import argparse
import bs4
import curses
import os
import re
import time

from contextlib import contextmanager
from datetime import datetime
from subprocess import check_output

@contextmanager
def newscreen():
    try:
        stdscr = curses.initscr()  # initialise it
        stdscr.clear()  # Clear the screen
        yield stdscr
    except KeyboardInterrupt:
        pass
    finally:
        curses.nocbreak()
        stdscr.keypad(0)
        curses.echo()
        curses.endwin()

def get_qstat(users, separate_arrays):
    x = check_output(['qstat', '-xt' if separate_arrays else '-x'])
    x = re.sub(b'\n', b'', x)
    soup = bs4.BeautifulSoup(x, 'lxml')

    if users:
        jobs = [x for x in soup.data.children if x.euser.string in users]
    else:
        jobs = soup.data.children

    out = []
    for job in jobs:
        if job.job_state.string == 'R':
            out.append(
                [
                job.job_id.string,
                job.job_name.string,
                job.euser.string,
                job.session_id.string if job.session_id is not None else'X',
                job.job_state.string,
                job.exec_host.string if job.exec_host is not None else'X',
                job.resources_used.walltime.string if job.resources_used is not None and job.resources_used.walltime is not None else'X',
                ])
        else:
            out.append(
                [
                job.job_id.string,
                job.job_name.string,
                job.euser.string,
                '-',
                job.job_state.string,
                '-',
                '-',
                ])
    return out


def print_qstat_to_screen(out, fstring1, fstring2, header, dstring):
    print(fstring1.format(*header))
    print(dstring)
    for _out in out:
        print(fstring2.format(*_out))

def print_qstat_to_curses(nsc, out, fstring1, fstring2, header, dstring, watch_duration):
    nsc.addstr(0, 0, datetime.now().strftime("%B %d, %Y  %H:%M:%S"), curses.A_REVERSE)
    nsc.addstr(1, 0, '')
    nsc.addstr(2, 0, 'qstat: refreshing Every {} seconds.'.format(watch_duration))
    nsc.addstr(3, 0, '')
    outstring = '\n'.join([fstring1.format(*header),
                           dstring,
                           '\n'.join([fstring2.format(*_out) for _out in out])
                           ])
    nsc.addstr(4, 0, outstring)
    nsc.refresh()


def qstat(users, watch_duration, separate_arrays):
    out = get_qstat(users, separate_arrays)

    header = ['JobID', 'JobName', 'Submitter', 'SessID', 'State', 'Cpu(s)', 'ElapsedTime']
    pads = [str(max([len(y[x]) for y in out + [header]])+2) for x in range(len(out[0]))]
    fstring1 = '{:^' + 's}|{:^'.join(pads) + 's}'
    fstring2 = '{:<' + 's}|{:<'.join(pads[:2]) + 's}|{:^' + 's}|{:^'.join(pads[2:]) + 's}'
    dstring = '-'*(sum([int(x) for x in pads]) + 20)  # (7) * 2 for `  ` + 6 * 1 for `|`

    if watch_duration == -1:
        print_qstat_to_screen(out, fstring1, fstring2, header, dstring)
    else:
        with newscreen() as nsc:
            while True:
                print_qstat_to_curses(nsc, out, fstring1, fstring2, header, dstring, watch_duration)
                time.sleep(watch_duration)
                out = get_qstat(users)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-U', '--user', dest='current_user', action='store_true',
                        help='Show only jobs submitted by current user.')
    parser.add_argument('-u', '--specific_users', dest='specific_users', type=str, nargs='+',
                        help='Show only jobs submitted by specific user(s). Negates -u.')
    parser.add_argument('-t', dest='t', action='store_true', help='Show individual array jobs.')
    parser.add_argument('-w', '--watch', dest='watch_duration', type=int, default=-1,
                        help='Refresh duration in seconds')
    params = parser.parse_args()

    assert params.watch_duration >= -1

    if params.specific_users:
        users = params.specific_users
    elif params.current_user:
        users = [os.environ["USER"]]
    else:
        users = []

    qstat(users, params.watch_duration, params.t)

if __name__ == '__main__':
    main()
