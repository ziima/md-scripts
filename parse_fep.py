"""Parse FEP output file."""
import os
import re
import shutil
import sys
from optparse import OptionParser

import VMD


PARSER = OptionParser(usage="vmd -e parse_fep.py -args forward backward")
PARSER.add_option("--output", default=".", help="Output directory [default: %default]")


def fix_fepout(filename, output):
    """
    Copies fepout to `output` directory and appends missing 'Free energy change' line
    """
    shutil.copy(filename, output)
    new_filename = os.path.join(output, os.path.basename(filename))

    last_step = None  # Last 'Free energy change' line
    fep_window = None  # Last 'NEW FEP WINDOW' line
    last_energy = None  # Last 'FepEnergy' line
    for line in open(new_filename):
        line = line.strip()
        if line.startswith('#Free energy change '):
            last_step = line
        elif line.startswith('#NEW FEP WINDOW: '):
            fep_window = line
        elif line.startswith('FepEnergy: '):
            last_energy = line
    assert last_step is not None
    assert fep_window is not None
    assert last_energy is not None

    # Find lambdas
    match = re.match(r'^#NEW FEP WINDOW: LAMBDA SET TO (?P<lambda>[-0-9.e]+) LAMBDA2 (?P<lambda2>[-0-9.e]+)$', fep_window)
    assert match is not None
    # XXX: Keep it as string. If converted to float, the fucking parsefep may fail.
    l, l2 = match.group('lambda'), match.group('lambda2')

    # Find delta
    chunks = last_energy.split()
    assert chunks[0] == 'FepEnergy:'
    assert chunks[-1]
    delta = float(chunks[-1])

    # Find net change
    match = re.search(r'net change until now is (?P<net_change>[-0-9.]+)$', last_step)
    assert match is not None
    net_change = float(match.group('net_change'))

    with open(new_filename, 'a') as out:
        out.write(
            '#Free energy change for lambda window [ %s %s ] is %.4f ; net change until now is %.4f\n'
            % (l, l2, delta, net_change + delta)
        )
    return new_filename


def main():
    options, args = PARSER.parse_args()

    if len(args) != 2:
        sys.stderr.write("This script requires exactly 2 arguments.")
        quit(1)

    forward, backward = args
    output = os.path.abspath(options.output)

    forward = os.path.abspath(forward)
    backward = os.path.abspath(backward)
    if not os.path.isdir(output):
        os.mkdir(output)
    os.chdir(output)

    forward = fix_fepout(forward, output)
    backward = fix_fepout(backward, output)

    VMD.evaltcl('package require parsefep')
    VMD.evaltcl('parsefep -forward %s -backward %s -bar' % (forward, backward))
    quit()


if __name__ == '__main__':
    main()
