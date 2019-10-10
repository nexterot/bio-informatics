import argparse

STATES = (DESC, DNA) = range(2)

DONE = 0
LEFT = 1
UP   = 2
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


def needleman_wunsch(a, b, g, scoring_func):
    D = [[0] * (len(b.dna) + 1) for _ in range(len(a.dna) + 1)]
    PTR = [[0] * (len(b.dna) + 1) for _ in range(len(a.dna) + 1)]
    D[0][0] = 1
    PTR[0][0] = DONE
    for i in range(1, len(a.dna) + 1):
        D[i][0] = i * g
        PTR[i][0] = UP
    for i in range(1, len(b.dna) + 1):
        D[0][i] = i * g
        PTR[0][i] = LEFT
    # dynamic programming step
    for i in range(1, len(a.dna) + 1):
        for j in range(1, len(b.dna) + 1):
            score_up = D[i-1][j] + g
            print(D[i-1][j], g)
            score_left = D[i][j-1] + g
            score_diag = D[i-1][j-1] + scoring_func(a.dna[i-1], a.dna[j-1])
            print(score_up, score_left, score_diag)
            score_max = max(score_up, score_left, score_diag)
            D[i][j] = score_max
            if score_max == score_up:
                PTR[i][j] = UP
            elif score_max == score_left:
                PTR[i][j] = LEFT
            else:
                PTR[i][j] = DIAG
    print("D=", D)
    print("PTR=", PTR)


def main():
    parser = argparse.ArgumentParser(description='NW')
    parser.add_argument('-i', metavar='file', type=str, nargs='+', help='input files')
    parser.add_argument('-o', type=str, nargs='+', help='output file')
    parser.add_argument('-g', type=int, nargs='+', help='gap')
    args = parser.parse_args()
    print(args)
    scoring_func = scoring_blosum62
    if len(args.i) == 1:
        r = parse_file(args.i[0], 2)
        print(r)
        needleman_wunsch(r[0], r[1], args.g[0], scoring_func)
    elif len(args.i) == 2:
        r = parse_file(args.i[0], 1) + parse_file(args.i[1], 1)
        needleman_wunsch(r[0], r[1], args.g[0], scoring_func)
        print(r)
    else:
        print("Provide 1 or 2 input files")
        return


if __name__ == "__main__":
    main()
