import argparse
import math

from sys import stdout

STATES = (DESC, DNA) = range(2)

DONE = 0
LEFT = 1
UP = 2
DIAG = 3

BLOSUM = (
    #  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
    (4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, -2, -1, 0, -4),  # A
    (-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, -1, 0, -1, -4),  # R
    (-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, 3, 0, -1, -4),  # N
    (-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, 4, 1, -1, -4),  # D
    (0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4),  # C
    (-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, 0, 3, -1, -4),  # Q
    (-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4),  # E
    (0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, -1, -2, -1, -4),  # G
    (-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, 0, 0, -1, -4),  # H
    (-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, -3, -3, -1, -4),  # I
    (-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, -4, -3, -1, -4),  # L
    (-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, 0, 1, -1, -4),  # K
    (-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, -3, -1, -1, -4),  # M
    (-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, -3, -3, -1, -4),  # F
    (-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, -2, -1, -2, -4),  # P
    (1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, 0, 0, 0, -4),  # S
    (0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, -1, -1, 0, -4),  # T
    (-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, -4, -3, -2, -4),  # W
    (-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, -3, -2, -1, -4),  # Y
    (0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, -3, -2, -1, -4),  # V
    (-2, -1, 3, 4, -3, 0, 1, -1, 0, -3, -4, 0, -3, -3, -2, 0, -1, -4, -3, -3, 4, 1, -1, -4),  # B
    (-1, 0, 0, 1, -3, 3, 4, -2, 0, -3, -3, 1, -1, -3, -1, 0, -1, -3, -2, -2, 1, 4, -1, -4),  # Z
    (0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0, 0, -2, -1, -1, -1, -1, -1, -4),  # X
    (-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, 1)  # *
)

DNA_FULL = (
    # A   T   G   C
    (5, -4, -4, -4),  # A
    (-4, 5, -4, -4),  # T
    (-4, -4, 5, -4),  # G
    (-4, -4, -4, 5),  # C
)


def get_blosum_index(c):
    if c == 'A':
        return 0
    if c == 'R':
        return 1
    if c == 'N':
        return 2
    if c == 'D':
        return 3
    if c == 'C':
        return 4
    if c == 'Q':
        return 5
    if c == 'E':
        return 6
    if c == 'G':
        return 7
    if c == 'H':
        return 8
    if c == 'I':
        return 9
    if c == 'L':
        return 10
    if c == 'K':
        return 11
    if c == 'M':
        return 12
    if c == 'F':
        return 13
    if c == 'P':
        return 14
    if c == 'S':
        return 15
    if c == 'T':
        return 16
    if c == 'W':
        return 17
    if c == 'Y':
        return 18
    if c == 'V':
        return 19
    if c == 'B':
        return 20
    if c == 'Z':
        return 21
    if c == 'X':
        return 22
    return 23


def get_dna_full_index(c):
    if c == 'A':
        return 0
    if c == 'T':
        return 1
    if c == 'G':
        return 2
    if c == 'C':
        return 3
    else:
        raise AssertionError


def scoring_dna_full(a, b):
    return DNA_FULL[get_dna_full_index(a)][get_dna_full_index(b)]


def scoring_blosum62(a, b):
    return BLOSUM[get_blosum_index(a)][get_blosum_index(b)]


def scoring_default(a, b):
    return 1 if a == b else -1


class Record:
    def __init__(self):
        self.desc = ""
        self.dna = ""

    def __repr__(self):
        return "Record(desc={},dna={})".format(self.desc, self.dna)


def parse_file(file_name, num_records):
    state = STATES[DESC]
    res = []
    with open(file_name) as f:
        rec = Record()
        for line in f.readlines():
            if num_records == 0:
                break
            if state == STATES[DESC]:
                if line.startswith(">"):
                    state = STATES[DNA]
                    rec.desc = line[1:].rstrip("\n")
            elif state == STATES[DNA]:
                if line.startswith(">"):
                    res.append(rec)
                    num_records -= 1
                    rec = Record()
                    rec.desc = line[1:].rstrip("\n")
                else:
                    line = line.rstrip("\n")
                    assert line.isalpha()
                    rec.dna += line
        if rec.desc and rec.dna:
            res.append(rec)
    return res


def pprint(s, m):
    print(s + "[")
    for row in m:
        print("\t", row, sep="")
    print("]")


def prints(d, s, file):
    stream = file if file else stdout
    print(d, file=stream)
    N = 100
    if len(s) <= N:
        print(s, file=stream)
        return
    num = len(s) // N
    for i in range(0, num):
        print(s[i:i+N], file=stream)
    print(s[num*N:], file=stream)


def smith_waterman(a, b, g, scoring_func, output_filename):
    i1 = 0
    j1 = 0
    max_val = -math.inf
    print(max_val)

    D = [[0] * (len(b.dna) + 1) for _ in range(len(a.dna) + 1)]
    PTR = [[0] * (len(b.dna) + 1) for _ in range(len(a.dna) + 1)]
    D[0][0] = 0
    PTR[0][0] = DONE
    for i in range(1, len(a.dna) + 1):
        D[i][0] = 0
        PTR[i][0] = UP
    for i in range(1, len(b.dna) + 1):
        D[0][i] = 0
        PTR[0][i] = LEFT
    # динамическое программирование
    for i in range(1, len(a.dna) + 1):
        for j in range(1, len(b.dna) + 1):
            score_up = D[i-1][j] + g
            score_left = D[i][j-1] + g
            score_diag = D[i-1][j-1] + scoring_func(a.dna[i-1], b.dna[j-1])
            # print(score_up, score_left, score_diag)
            score_max = max(0, score_up, score_left, score_diag)
            D[i][j] = score_max
            if score_max == score_up:
                PTR[i][j] = UP
            elif score_max == score_left:
                PTR[i][j] = LEFT
            else:
                PTR[i][j] = DIAG

            if score_max > max_val:
                print(score_max)
                max_val = score_max
                i1 = i
                j1 = j
    # pprint("D", D)
    # pprint("PTR", PTR)
    # обратный ход
    count = 0
    s = len(a.dna) + len(b.dna)
    res_a = ['?'] * s
    res_b = ['?'] * s

    i = i1
    j = j1

    while i != 0 and j != 0 and D[i][j] != 0:
        if PTR[i][j] == DIAG:
            res_a[count] = a.dna[i-1]
            res_b[count] = b.dna[j-1]
            count += 1
            i -= 1
            j -= 1
        elif PTR[i][j] == LEFT:
            res_a[count] = '-'
            res_b[count] = b.dna[j-1]
            count += 1
            j -= 1
        elif PTR[i][j] == UP:
            res_a[count] = a.dna[i-1]
            res_b[count] = '-'
            count += 1
            i -= 1
    for k in range(0, count // 2):
        tmp = res_a[k]
        res_a[k] = res_a[count - 1 - k]
        res_a[count - 1 - k] = tmp
        tmp = res_b[k]
        res_b[k] = res_b[count - 1 - k]
        res_b[count - 1 - k] = tmp
    score = D[i1][j1]
    a1 = "".join(res_a).rstrip("?")
    a2 = "".join(res_b).rstrip("?")
    if output_filename:
        with open(output_filename, "w") as f:
            prints("seq1: ", a1, f)
            prints("seq2: ", a2, f)
            print("score:", score, file=f)
        return
    prints("seq1: ", a1, None)
    prints("seq2: ", a2, None)
    print("score:", score)


def main():
    parser = argparse.ArgumentParser(description='NW')
    parser.add_argument('-i', metavar='file', type=str, nargs='+', help='input files')
    parser.add_argument('-o', type=str, help='output file')
    parser.add_argument('-s', type=str, help='')
    parser.add_argument('-g', type=int, help='gap')
    args = parser.parse_args()
    if args.g is None:
        args.g = -2
    # print(args)
    scoring_func = None
    if args.s == "blosum62":
        scoring_func = scoring_blosum62
    elif args.s == "dna_full":
        scoring_func = scoring_dna_full
    else:
        scoring_func = scoring_default

    if len(args.i) == 1:
        r = parse_file(args.i[0], 2)
        # print(r)
        smith_waterman(r[0], r[1], args.g, scoring_func, args.o)
    elif len(args.i) == 2:
        r = parse_file(args.i[0], 1) + parse_file(args.i[1], 1)
        # print(r)
        smith_waterman(r[0], r[1], args.g, scoring_func, args.o)
    else:
        print("invalid number of input files")
        return


if __name__ == "__main__":
    main()
