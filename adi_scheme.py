def adi_scheme(initial_temp, r_temp, l_temp, top_temp, bottom_temp, time_step, step_size, length, time, diffusivity):
    """
    The alternating-direction implicit, or ADI, scheme provides a means for solving parabolic
    equations in two spatial dimensions using tri-diagonal matrices.
    The following inputs are an example of a specific case. Taken from Steven C. Chapra, Numerical methods for
    engineers, Finite difference: Parabolic equations.
    in_t = [[75, 75, 75, 75, 75], [0, 0, 0, 0, 100], [0, 0, 0, 0, 100], [0, 0, 0, 0, 100], [50, 50, 50, 50, 50]]
    in_t is the initial temperature distribution, which includes the boundary temperatures.
    time_step = 10 sec
    step_size = 10 cm
    length = 40 cm
    time = 300 sec
    diffusivity = 0.835 cm^2/s.#  diffusivity is the thermal diffusivity of the material.
    """
    # To create an array in_t with  initial conditions  and the boundary temperatures
    x_n = int(length/step_size) if step_size != 0 else 0
    left_side_temp = [l_temp for _ in range(x_n + 1)]
    right_side_temp = [r_temp for _ in range(x_n + 1)]
    in_t = []  # to initialize the array
    for i in range(x_n+1):
        if i == 0:
            in_t.append(left_side_temp)
        elif i == x_n:
            in_t.append(right_side_temp)
        else:
            temp = []  # temporary array to hold the top temp, interior nodes temp, and bottom temperatures
            for j in range(x_n + 1):
                if j == 0:
                    temp.append(bottom_temp)
                elif j == x_n:
                    temp.append(top_temp)
                else:
                    temp.append(initial_temp)
            in_t.append(temp)
    la = diffusivity * time_step / (step_size ** 2) if step_size != 0 else 0  # constant lamda value
    # coefficient matrix a
    a = [[0 for _ in range(1, len(in_t) - 1)] for _ in range(1, len(in_t) - 1)]
    for i in range(1, len(in_t) - 1):
        for j in range(1, len(in_t) - 1):
            if j == i:
                a[i - 1][j - 1] = 2 * (1 + la)
            elif j == i - 1 or j == i + 1:
                a[i - 1][j - 1] = -la
    # decompose matrix a into lower and upper triangle matrices
    for i in range(1, len(a)):
        for j in range(1, len(a)):
            if j == i:
                a[i][j - 1] = a[i][j - 1] / a[i - 1][j - 1] if a[i - 1][j - 1] != 0 else 0
                a[i][j] = a[i][j] - a[i][j - 1] * a[i - 1][j]
    # for iteration over several time step
    # to store over each time iteration initialize temp_over_time
    temp_over_time = []
    for t in range(int(time / time_step)):
        # for t = l to t = l+1/2 y_implicit
        b = [0 for _ in range(len(a))]
        for r in range(1, len(in_t) - 1):  # remember to add additional loop for T[(2,js)] and T[(3,js)]
            for k in range(1, len(in_t) - 1):
                if k == 1:
                    b[k - 1] = la * in_t[r - 1][k] + 2 * (1 - la) * in_t[r][k] + la * in_t[r + 1][k] + la * in_t[r][
                        k - 1]
                elif k == len(in_t) - 2:
                    b[k - 1] = la * in_t[r - 1][k] + 2 * (1 - la) * in_t[r][k] + la * in_t[r + 1][k] + la * in_t[r][
                        k + 1]
                else:
                    b[k - 1] = la * in_t[r - 1][k] + 2 * (1 - la) * in_t[r][k] + la * in_t[r + 1][k]
            # Forward substitution to modify b
            for i in range(1, len(b)):
                b[i] = b[i] - b[i - 1] * a[i][i - 1]
            # Backward substitution to find T(1,js):
            temp_col = [0 for _ in range(len(b))]
            temp_col[-1] = b[-1] / a[-1][-1] if a[-1][-1] != 0 else 0
            for i in range(len(temp_col) - 2, -1, -1):
                temp_col[i] = (b[i] - a[i][i + 1] * temp_col[i + 1]) / a[i][i] if a[i][i] != 0 else 0
            # update in_t[1][j's]
            for j in range(1, len(a) + 1):
                in_t[r][j] = temp_col[j - 1]
        # for t = l+1/2 to t= l+1 x_implicit
        b = [0 for _ in range(len(a))]
        temp_instant = []  # to hold a temperature for one time_step iteration and append to temp_over_time
        for r in range(1, len(in_t) - 1):  # remember to add additional loop for T[(is,2)] and T[(is,3)]
            for k in range(1, len(in_t) - 1):
                if k == 1:
                    b[k - 1] = la * in_t[k][r - 1] + 2 * (1 - la) * in_t[k][r] + la * in_t[k][r + 1] + la * in_t[k - 1][
                        r]
                elif k == len(in_t) - 2:
                    b[k - 1] = la * in_t[k][r - 1] + 2 * (1 - la) * in_t[k][r] + la * in_t[k][r + 1] + la * in_t[k + 1][
                        r]
                else:
                    b[k - 1] = la * in_t[k][r - 1] + 2 * (1 - la) * in_t[k][r] + la * in_t[k][r + 1]
            # Forward substitution to modify b
            for i in range(1, len(b)):
                b[i] = b[i] - b[i - 1] * a[i][i - 1]
            # Backward substitution to find the temp row i.e, T(1,1),T(2,1),T(3,1) ...
            temp_row = [0 for _ in range(len(b))]
            temp_row[-1] = b[-1] / a[-1][-1] if a[-1][-1] != 0 else 0
            for i in range(len(temp_row) - 2, -1, -1):
                temp_row[i] = (b[i] - temp_row[i + 1] * a[i][i + 1]) / a[i][i] if a[i][i] != 0 else 0
            # update in_t[i's][j]
            for i in range(1, len(a) + 1):
                in_t[i][r] = temp_row[i - 1]
            temp_instant.append(temp_row)
        temp_over_time.append(temp_instant)
    return temp_over_time


case = adi_scheme(0, 50, 75, 100, 0, 10, 10, 40, 300, 0.835)
print(case)
