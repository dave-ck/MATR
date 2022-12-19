# MATR script wirtten by David Kutner, commissioned by Richard Harris, both of Durham University. Free to use with attribution.
import pandas as pd
import numpy as np
from scipy.optimize import linprog


def slackify_guess(A_eq, b_eq, x0):
    # we create a new "slack" variable si for each variable xi; we denote the guess constant x0[i] by gi in comments
    nEqs = len(b_eq)
    nVars = A_eq.shape[1]
    assert (A_eq.shape[0] == nEqs)
    # extend A_eq to include (with factor zero) the slack variables in the equalities
    A_eq = np.hstack((A_eq, np.zeros((nEqs, nVars))))
    # inequalities (greather-than->gt) to set xi-si<=gi and -xi-si<=-gi (i.e. if we minimize si it is |xi-gi|)
    A_gt = np.zeros((nVars * 2, nVars * 2))
    current_eq = 0
    for varIndex in range(nVars):
        A_gt[current_eq, varIndex] = 1
        A_gt[current_eq, varIndex + nVars] = -1
        current_eq += 1
        A_gt[current_eq, varIndex] = -1
        A_gt[current_eq, varIndex + nVars] = -1
        current_eq += 1
    # duplicate x0, so b=(g0,g0,g1,g1,...) then alternatively negate so b=(g0,-g0,g1,-g1,...)
    b_gt = np.repeat(x0, 2) * np.tile([1, -1], nVars)
    c = np.repeat(np.array([0, 1], dtype=float), nVars)
    print(nEqs, nVars)
    return c, A_eq, b_eq, A_gt, b_gt


def generate_matr_problem(matrix_shape, suppress_prob):
    """
    :rtype: object
    :param matrix_shape: desired matrix shape e.g. (3,7)
    :param suppress_prob: probability of suppressed values e.g. 0.15
    :return: ground truth matrix and input matrix, as a tuple of 2 np.ndarrays
    """
    # populate a matrix with uniform random entries in [1,1000]
    ground_truth = np.random.randint(1, 1000, matrix_shape).astype(dtype=np.float64)
    # calculate column, row and grand totals
    ground_truth = np.append(ground_truth, ground_truth.sum(axis=0).reshape(1, matrix_shape[1]), axis=0)
    ground_truth = np.append(ground_truth, ground_truth.sum(axis=1).reshape(matrix_shape[0] + 1, 1), axis=1)

    n, m = matrix_shape
    matrix_shape = n + 1, m + 1
    # create boolean "mask" matrix to suppress entries of input
    suppress = np.random.random(matrix_shape) < suppress_prob
    suppress[-1, :] = False
    suppress[:, -1] = False
    input_matrix = ground_truth.copy()
    # set suppressed cells to nan (as when a blank cell is read with pandas)
    input_matrix[suppress] = np.nan

    guess_error = 0.5 + np.random.random(matrix_shape)  # guess is between 50% and 150% of real value
    guess_matrix = np.multiply(ground_truth, guess_error)
    return ground_truth, input_matrix, guess_matrix


def solve_matr_problem(input_npa, guess_npa=None):
    # separate out the row and column totals (ignore grand total)
    column_totals = input_npa[-1, :-1]
    row_totals = input_npa[:-1, -1]

    # From here we refer to:
    #   the m*n input and guess matrices,
    #   the m*1 column total, and
    #   the 1*n row totals
    m = input_npa.shape[0] - 1
    n = input_npa.shape[1] - 1

    # iterate over entire array, obtain LP_vars
    blank_coords = list()
    for x in range(n):
        for y in range(m):
            # print(x, y, "in loop - testing - v[y,x]=", input_npa[y, x])
            if np.isnan(input_npa[y, x]):
                blank_coords.append((y, x))

    # initialize matrix A_eq and associated variables
    N_blanks = len(blank_coords)
    N_vars = N_blanks  # for each blank (at y,x), we have 3 variables: v_(y,x), plus_(y,x), and minus(y_x)

    N_eqs = m + n
    current_eq = 0
    A_eq = np.zeros((N_eqs, N_vars))
    b_eq = np.zeros(N_eqs)

    # guess vector x0
    x0 = np.zeros(len(blank_coords))
    for i in range(len(blank_coords)):
        coord = blank_coords[i]
        x0[i] = guess_npa[coord]

    # iterate over rows, find unallocated amount
    for y in range(m):
        row_unallocated = row_totals[y]
        for x in range(n):
            if np.isnan(input_npa[y, x]):
                current_var = blank_coords.index((y, x))
                A_eq[current_eq, current_var] = 1  # v_(y,x) is a variable in our LP
            else:
                row_unallocated -= input_npa[y, x]
        b_eq[current_eq] = row_unallocated  # sum of our variables should equal this unallocated amout
        current_eq += 1  # increment current equation

    # iterate over cols, find unallocated amount; like above
    for x in range(n):
        column_unallocated = column_totals[x]
        for y in range(m):
            if np.isnan(input_npa[y, x]):
                current_var = blank_coords.index((y, x))
                A_eq[current_eq, current_var] = 1
            else:
                column_unallocated -= input_npa[y, x]
        b_eq[current_eq] = column_unallocated
        current_eq += 1

    # Solve the linear program using https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linprog
    slackify_guess(A_eq, b_eq, x0)
    c, A_eq, b_eq, A_gt, b_gt = slackify_guess(A_eq, b_eq, x0)
    sol = linprog(c, A_eq=A_eq, b_eq=b_eq, A_ub=A_gt, b_ub=b_gt)
    x = sol.x[:N_blanks]

    output_npa = input_npa.copy()

    for i in range(N_blanks):
        blank_coord = blank_coords[i]
        value_for_blank = x[i]
        assert (np.isnan(output_npa[blank_coord]))  # ensure np.nan only value overwritten
        output_npa[blank_coord] = value_for_blank
    output_npa[-1, :] = output_npa[:-1, :].sum(axis=0)  # update sums
    output_npa[:, -1] = output_npa[:, :-1].sum(axis=1)
    out_column_totals = input_npa[-1, :-1]
    out_row_totals = input_npa[:-1, -1]
    column_error = np.abs(column_totals - out_column_totals).sum()
    row_error = np.abs(row_totals - out_row_totals).sum()
    distance_from_guess = np.abs(x0 - x).sum()
    total_error = column_error + row_error
    print("Total error:", total_error, "Distance from guess:", distance_from_guess)
    return output_npa


if __name__ == "__main__":
    input_df = pd.read_excel("input.xlsx", header=None)
    guess_df = pd.read_excel("guess.xlsx", header=None)
    input_npa = input_df.to_numpy()
    guess_npa = guess_df.to_numpy()
    output_npa = solve_matr_problem(input_npa, guess_npa)
    output_df = pd.DataFrame(output_npa)
    output_df.to_excel("output.xlsx", header=None, index=False)
