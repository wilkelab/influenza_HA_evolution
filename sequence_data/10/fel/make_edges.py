import sys, os

def main():
    in_filename = sys.argv[1]
    out_filename = sys.argv[2]
    make_edges(in_filename, out_filename)

def make_edges(in_filename, out_filename):
    positive = []
    lines = list(open(in_filename, 'r').readlines())
    for line in enumerate(lines):
        if line[0] == 0:
            continue
        else:
            split_line = line[1].split(',')
            if float(split_line[0]) > 1 and float(split_line[2]) < 0.01 and (int(line[0]) - 16) > 0:
                positive.append(int(line[0]) - 16)

    print(positive)
    out_file = open(out_filename, 'w')
    #out_file.write('n1,n2\n')
    for site1 in range(len(positive) - 1):
        for site2 in range(site1 + 1, len(positive)):
            out_file.write(str(positive[site1]) + ',' + str(positive[site2]) + '\n')
    out_file.close()

main()
