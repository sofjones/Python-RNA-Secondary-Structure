'''
Pairs map to RNA pairs 'G' 'C' strongest pair
A replaced with T because using a DNA dataset
'''
score_1 = [['G', 'C'], ['C', 'G']]
score_2 = [['A', 'T'], ['T', 'A']]
score_3 = [['G', 'T'], ['T', 'G']]


# initialize square matrix
def initMatrix(size):
    mat = [[0 for row in range(size)] for col in range(size)]
    return mat

#score pair make sure pairs are at least the minimum loop apart
def scorePair(i, j, seq, min_loop):
    pair = [seq[i], seq[j]]
    bonus = 0
    if i >= j - (min_loop * 1.5):
        bonus = -3
    elif i >= j - min_loop:
        return 0
    if (pair in score_1):
        score = 2 + bonus
    elif (pair in score_2):
        score = 1 + bonus
    elif (pair in score_3):
        score = 0 + bonus
    else:
        score = 0
    return score


#split up to avoid complexity in findBestCase
def case_four(i, j, mat):
    temp = 0
    for idx in range(i + 1, j - 1):
        temp = max(temp, (mat[i][idx] + mat[idx + 1][j]))
    return temp

#finds the best case out of 4 nussinov cases
def findBestCase(i, j, mat, seq, min_loop, energy_min=True):
    case1 = mat[i + 1][j]
    case2 = mat[i][j - 1]
    case3 = mat[i + 1][j - 1]
    if energy_min:
        case3 += scorePair(i, j, seq, min_loop)
    else:
        case3 += nuss_score(i, j, seq)
    case4 = case_four(i, j, mat)
    max_case = max(case1, case2, case3, case4)
    return max_case


#traceback runs until i is equal to j, or half the matrix is run through
def traceback(i, j, mat, seq, struc, min_loop, energy_min=True):
    if i >= j:
        return
    if energy_min == False:
        case3 = mat[i + 1][j - 1] + nuss_score(i, j, seq)
    else:
        case3 = mat[i + 1][j - 1] + scorePair(i, j, seq, min_loop)
        if addStacked(i, j, mat, seq, min_loop):
            case3 += 2
    if mat[i][j] == mat[i + 1][j]:
        traceback(i + 1, j, mat, seq, struc, min_loop, energy_min)
    elif mat[i][j] == mat[i][j - 1]:
        traceback(i, j - 1, mat, seq, struc, min_loop, energy_min)
    elif mat[i][j] == case3 or mat[i][j] == case3 + 2:
        struc.append((i + 1, j + 1))
        traceback(i + 1, j - 1, mat, seq, struc, min_loop, energy_min)
    else:
        for idx in range(i + 1, j - 1):
            if mat[i][j] == mat[i][idx] + mat[idx + 1][j]:
                traceback(i, idx, mat, seq, struc, min_loop, energy_min)
                traceback(idx + 1, j, mat, seq, struc, min_loop, energy_min)
                break

#creates dot bracket structure from base pairs
def dotBracket(seq, base_pairs, max_score):
    d_b = seq
    for i in range(0, len(seq)):
        open_par = [(x, y) for x, y in base_pairs if x == i + 1]
        closed_par = [(x, y) for x, y in base_pairs if y == i + 1]
        if open_par:
            d_b[i] = "("
        elif closed_par:
            d_b[i] = ")"
        else:
            d_b[i] = "."

    return d_b

# method to check if pairs are stacked
def addStacked(i, j, mat, seq, min_loop):
    past_score = scorePair(i, j, seq, min_loop)
    score = mat[i][j]
    if i >= j - (min_loop * 1.5):
        return False
    if past_score > 1 and i > 0 and mat[i + 1][j - 1] > 0:
        score += 2
        return True
    return False

#ensures that there is a one to one map from i to j
def check_results(results):
    open_par = 0
    closed_par = 0
    for c in results:
        if c == '(':
            open_par += 1
        elif c == ')':
            closed_par += 1
    if open_par - closed_par == 0:
        return True
    else:
        return False

def nuss_energy_min(seq, min_loop):
    length = len(seq)
    matrix = initMatrix(length)
    base_pairs = []
    # fill matrix energy min
    for idx in range(1, length):
        # only half matrix needs to be completed
        for i in range(length - idx):
            j = idx + i
            matrix[i][j] = findBestCase(i, j, matrix, seq, min_loop)

    for idx in range(1, length):
        for i in range(1, length - idx):
            j = idx + i
            is_stacked = addStacked(i, j, matrix, seq, min_loop)
            if is_stacked:
                base_pairs.append((i, j))
    traceback(i, j, matrix, seq, base_pairs, min_loop)
    # max score shown in top right hand corner of matrix
    max_score = matrix[0][length - 1]
    results = dotBracket(seq, base_pairs, max_score)
    correct_output = check_results(results)
    while (correct_output == False):
        base_pairs.pop(0)
        results = dotBracket(seq, base_pairs, max_score)
        correct_output = check_results(results)
    return results

def energy_min(seq):
    min_loop = len(seq) * 0.25
    nuss_MFE = nuss_energy_min(seq, min_loop)
    return nuss_MFE
