'''
included for reference 
this is an early version of our Nussinov algorithm
which has since been replaced by the file energy_min.py
''''

import sys
import matplotlib.pyplot as plt
import forgi.visual.mplotlib as fvm
import forgi


def build_matrix(sequence):
    size = len(sequence)
    matrix = []
    matrix.append([" "] + sequence)
    for i in range(size):
        matrix.append([sequence[i]])
        for j in range(size):
            matrix[i+1].append(0)
    return matrix

def complement(n1, n2):
    if n1 == "a" and n2 in "ut":
        return 1
    elif  n1 in "tu" and n2 == "a":
        return 1
    elif  n1 in "ut" and n2 == "g":
        return 1
    elif  n1 == "g" and n2 in "ut":
        return 1
    elif  n1 == "c" and n2 == "g":
        return 1
    else:
        return 1 if (n1 == "g" and n2 == "c") else 0

def fill_matrix(matrix):
    counter = 1
    for i in range(1, len(matrix)):
        for j in range(i + 1, len(matrix)):
            matrix[j-i][j] = nuss_func(matrix, j-i, j)


def traceback(results, matrix, i, j):
    # print("in", i, j)
    if j <= i:
        return
    elif matrix[i][j] == matrix[i][j-1]:
        traceback(results, matrix, i, j-1)
        return
    else:
        for k in range(i, j):
            if complement(matrix[0][k], matrix[0][j]):
                if k > 1:
                    if matrix[i][j] == matrix[i][k-1] + matrix[k+1][j-1] + 1:
                        results.append((k, j))
                        traceback(results, matrix, i, k-1)
                        traceback(results, matrix, k+1, j-1)
                        return
                        # print("out3")
                elif matrix[i][j] == matrix[k+1][j-1] + 1:
                    results.append((k, j))
                    traceback(results, matrix, i, k-1)
                    traceback(results, matrix, k+1, j-1)
                    return
                    # print("out3")


def nuss_func(matrix, i, j):
    # hanging_i = matrix[1+1][j]
    hanging_j = matrix[i][j-1]
    # joined = matrix[i+1][j-1] + complement(matrix[0][j], matrix[i][0])
    sub_structures = get_max_substructures(matrix, i, j)
    # return max(hanging_i, hanging_j, joined, sub_structures)
    return max(hanging_j, sub_structures)


def get_max_substructures(matrix, i, j):
    max = 0
    temp = 0
    for k in range(i, j):
        if complement(matrix[0][k], matrix[0][j]):
            if k > 1 :
                temp = matrix[i][k-1] + matrix[k+1][j-1] + 1
            else:
                temp = matrix[k+1][j-1] + 1
        if temp > max:
            max = temp
    return max


def nussinov(sequence):
    # print(sequence)
    matrix = build_matrix(sequence)
    fill_matrix(matrix)
    results = []
    traceback(results, matrix, 1, len(matrix[0])-1)
    return results


def print_matrix(matrix):
    print('\n'.join([''.join(['%4s' % item for item in row])for row in matrix]))

