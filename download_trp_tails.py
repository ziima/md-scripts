"""Download tail segments of all human TRP channels."""
import os

from urllib import quote_plus
from urllib2 import urlopen

from xml.etree import ElementTree


OUTPUT = os.path.join(os.path.dirname(__file__), 'trp_tails.fasta')

DIRECT_URL = "http://www.uniprot.org/uniprot/%(code)s.xml"
QUERY_URL = 'https://www.uniprot.org/uniprot/?query=%(query)s&format=xml'
QUERY = 'name:"Transient receptor potential"' \
        ' AND reviewed:yes' \
        ' AND organism:"Homo sapiens"'
#         ' AND family:"1.A.4"' \
XML_NS = '{http://uniprot.org/uniprot}'


# How many residues around transmembrane region should be left
BEGIN_MARGIN = 5
END_MARGIN = 40


def get_code(entry):
    """
    Returns entry accession code.
    """
    return entry.find('%saccession' % XML_NS).text


def process_entry(entry):
    """Process <entry> element.

    Return pair name-sequence or None.
    """
    name = entry.find('%sname' % XML_NS).text
    print "Name:", name

    organism = entry.find("%(ns)sorganism/%(ns)sname[@type='scientific']" % {'ns': XML_NS}).text
    print "Organism:", organism

    sequence = entry.find('%ssequence' % XML_NS).text.replace('\n', '')

    features = entry.findall('%sfeature' % XML_NS)
    transmembrane_regions = [f for f in features if f.attrib['type'] == 'transmembrane region']
    print "Found %d transmembrane regions" % len(transmembrane_regions)

    if not transmembrane_regions:
        print "Skipping"
        return None

    # Find beginning and end of transmembrane region
    begin = int(transmembrane_regions[-1].find('%(ns)slocation/%(ns)sbegin' % {'ns': XML_NS}).attrib['position'])
    end = int(transmembrane_regions[-1].find('%(ns)slocation/%(ns)send' % {'ns': XML_NS}).attrib['position'])
    if name == 'TRPM1_HUMAN':
        # TRPM1 has some additional transmembrane region and the end
        begin = int(transmembrane_regions[-2].find('%(ns)slocation/%(ns)sbegin' % {'ns': XML_NS}).attrib['position'])
        end = int(transmembrane_regions[-2].find('%(ns)slocation/%(ns)send' % {'ns': XML_NS}).attrib['position'])
    print "Range from %d to %d (1-based)" % (begin, end)

    # Compute slice boundaries
    begin = max(begin - BEGIN_MARGIN - 1, 0)
    end = min(len(sequence), end + END_MARGIN)
    print "Slice: %d:%d (0-based)" % (begin, end)

    return name, sequence[begin:end]


results = {}

print "Query:", QUERY
url = QUERY_URL % {'query': quote_plus(QUERY)}
print url
response = urlopen(url)
tree = ElementTree.parse(response)
for entry in tree.findall('%sentry' % XML_NS):
    code = get_code(entry)
    print "=" * 10, code, "=" * 10
    if code in results:
        print "SKIP: Already processed"
        continue
    result = process_entry(entry)
    if result is not None:
        results[code] = result


CUSTOM_SEQUENCES = (
    'A8EVM5',  # Bacterial Voltage gated channel (pdb: 4EKW)
    'A0L5S6',  # Bacterial Voltage gated channel (pdb: 4F4L)
    'Q9P0L9',  # PKD2L1/TRPP3
    'Q9NZM6',  # PKD2L2/TRPP5
)
for code in CUSTOM_SEQUENCES:
    response = urlopen(DIRECT_URL % {'code': code})
    tree = ElementTree.parse(response)
    for entry in tree.findall('%sentry' % XML_NS):
        code = get_code(entry)
        print "=" * 10, code, "=" * 10
        if code in results:
            print "SKIP: Already processed"
            continue
        result = process_entry(entry)
        if result is not None:
            results[code] = result


def sort_by_name(item):
    return item[1][0]


print "=" * 30
print "Writing output"
with open(OUTPUT, 'w') as output:
    for (code, (name, sequence)) in sorted(results.items(), key=sort_by_name):
        # Append sequence to output
        output.write('>%(name)s\n' % {'name': name})
        output.write(sequence)
        output.write('\n')

print "Found %d sequences" % len(results)
