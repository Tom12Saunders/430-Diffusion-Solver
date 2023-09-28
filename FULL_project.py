import numpy as np
import matplotlib.pyplot as plt
import sys

##########################################  PROBLEM 1  ##########################################################################

class UTF8Encoder:
    def __init__(self, std):
        self.std = std

    def write(self, message):
        self.std.buffer.write(message.encode('utf-8'))

    def flush(self):
        pass  # Depending on use, you may want to implement

original_stdout = UTF8Encoder(sys.stdout)
output = 'Problem_1_output.txt'
with open (output, 'w') as theSack:
    # Redirecting the standard output to the file
    original_stdout = sys.stdout  # Save the original standard output
    sys.stdout = UTF8Encoder(theSack)  # Redirect the output to the file using the UTF8Encoder

    # Constants and Input Parameters
    Sigma_t = 0.83  # cm^-1
    Sigma_a = 0.1  # cm^-1
    vSigma_f = 0.125  # cm^-1
    xB = 5.5  # cm
    n = 40
    h = xB / n
    tolerance = 1e-6
    D = 1 / (3 * Sigma_t)

    x = np.linspace(0, xB, n)  # Adjusted linspace to go from 0 to xB

    # Neutron speed
    v = 1.0  # cm/s

    # Analytical solution
    k_analytical = (4 * xB ** 2 * vSigma_f) / (D * np.pi ** 2 + 4 * Sigma_a * xB ** 2)
    a = np.sqrt((vSigma_f - k_analytical * Sigma_a) / (k_analytical * D))
    phi_analytical = np.cos(a * x)

    # Matrix K (Diffusion and Absorption)
    K = np.diag((2 * D / h ** 2 + Sigma_a) * np.ones(n))
    K += np.diag((-D / h ** 2) * np.ones(n - 1), -1)
    K += np.diag((-D / h ** 2) * np.ones(n - 1), 1)
    F = np.diag(vSigma_f * np.ones(n))

    # Reflective boundary at x = 0
    K[0, 0] = 2 * D / h + Sigma_a  # Changed from 2 * D / h ** 2 + Sigma_a
    K[0, 1] = -2 * D / h  # Changed from -D / h ** 2

    # Zero flux and zero value boundary at x = xB
    K[-1, :] = 0  # Resetting the last row of K
    K[-1, -1] = 1  # Setting the diagonal entry to 1 for the last row

    K_inv = np.linalg.inv(K)

    # Analytical a-value for characteristic shape
    a = np.sqrt((vSigma_f - k_analytical * Sigma_a) / (k_analytical * D))

    # Power Iteration
    phi = np.cos(a*x)  # Starting guess
    phi /= np.sum(phi * h)  # Normalize initial flux to represent a single neutron in the system

    k_vals = []  # To keep track of k-values through iterations
    L2_diffs = []  # To store the L2 norms of the differences

    while True:
        # Compute fission source term
        F_phi = vSigma_f * phi

        # Get phi_new using K_inv and current fission source
        phi_new = np.linalg.solve(K, F_phi)

        # Estimate k
        k_new = np.sum(phi_new) / np.sum(phi)
        k_vals.append(k_new)

        # Normalize phi_new
        phi_new /= np.sum(phi_new * h)

        # Compute the L2 norm of the difference
        L2_diff = np.linalg.norm(phi_new - phi, 2)
        L2_diffs.append(L2_diff)

        # Check for convergence
        if len(L2_diffs) > 1 and L2_diffs[-1] < tolerance and abs(1 / k_vals[-1] * (k_vals[-1] - k_vals[-2])) < tolerance:
            break

        # Update phi for next iteration
        phi = phi_new

    # Output
    print(f'\n********************Computation of k-value and Flux for Problem 1********************\n')
    print('[Number of Iterations]:', len(k_vals), f'\n[Number of Discrete Spatial Bins]: {n}\n\n[Boundary Conditions]: \n(Left) Reflective\n(Right) Zero-Flux\n')
    print(f'[Numerical k-value after {len(k_vals)} Iterations]: {k_vals[-1]:.5f}')
    print(f'[Analytical k-value]: {k_analytical:.5f}\n')
    print(f'[Difference between Analytical and Numerical k-values]:\n', k_vals[-1]-k_analytical)
    print('\n[Numerical k-values for Each Iteration]:', k_vals)
    print(f'\n[Numerical flux φ [neutrons⋅s⁻¹⋅cm⁻²] at (x = 0)]: {phi[0]:.5f}')
    print(f'[Numerical flux φ [neutrons⋅s⁻¹⋅cm⁻²] at (x=5.5)]: {phi[-1]:.5f}')
    print(f'\n[Final φ {n}-vector]:\n', phi)
    print('\n[The sum of neutrons from x=0 to x=5.5 (cm) for all times in the system]:', np.sum(phi * h))

# Completing printout
sys.stdout = original_stdout  # Resetting the standard output back to the console
print(f'Output written to {output}')

# Plotting if packages and necessary software installed
plt.plot(x, phi, label='Numerical')
plt.plot(x, phi_analytical, '--', label='Analytical')
plt.xlabel('x (cm)')
plt.ylabel('φ(x) [neutrons⋅s⁻¹⋅cm⁻²]')
plt.legend()
plt.title('Eigenfunctions')
plt.show()

##########################################  PROBLEM 2  ##########################################################################

class UTF8Encoder:
    def __init__(self, std):
        self.std = std

    def write(self, message):
        self.std.buffer.write(message.encode('utf-8'))

    def flush(self):
        pass  # Depending on use, you may want to implement

original_stdout = UTF8Encoder(sys.stdout)

with open('Problem_2_output.txt', 'w', encoding='utf-8') as theSack:
    sys.stdout = UTF8Encoder(theSack)  # Redirect the output to the file using the UTF8Encoder
    # Imaginatively-named function
    def compute_flux_eigenvalue(N):
        # Constants and Input Parameters
        Sigma_t = 0.83  # cm^-1
        Sigma_a = 0.1  # cm^-1
        vSigma_f = 0.125  # cm^-1
        xB = 5.5  # cm
        n = N
        h = xB / n
        tolerance = 1e-8
        D = 1 / (3 * Sigma_t)

        x = np.linspace(0, xB, n)  # Adjusted linspace to go from 0 to xB

        # Neutron speed
        v = 1.0  # cm/s

        # Analytical solution
        k_analytical = (4 * xB ** 2 * vSigma_f) / (D * np.pi ** 2 + 4 * Sigma_a * xB ** 2)
        a = np.sqrt((vSigma_f - k_analytical * Sigma_a) / (k_analytical * D))
        phi_analytical = np.cos(a * x)

        # Matrix K (Diffusion and Absorption)
        K = np.diag((2 * D / h ** 2 + Sigma_a) * np.ones(n))
        K += np.diag((-D / h ** 2) * np.ones(n - 1), -1)
        K += np.diag((-D / h ** 2) * np.ones(n - 1), 1)
        F = np.diag(vSigma_f * np.ones(n))

        # Reflective boundary at x = 0
        K[0, 0] = 2 * D / h + Sigma_a
        K[0, 1] = -2 * D / h

        # Zero flux and zero value boundary at x = xB
        K[-1, :] = 0  # Resetting the last row of K
        K[-1, -1] = 1  # Setting the diagonal entry to 1 for the last row

        K_inv = np.linalg.inv(K)

        # Analytical a-value for characteristic shape
        a = np.sqrt((vSigma_f - k_analytical * Sigma_a) / (k_analytical * D))

        # Power Iteration
        phi = np.cos(a * x)  # Starting guess
        phi /= np.sum(phi * h)  # Normalize initial flux to represent a single neutron in the system

        k_vals = []  # To keep track of k-values through iterations
        L2_diffs = []  # To store the L2 norms of the differences

        while True:
            # Compute fission source term
            F_phi = vSigma_f * phi

            # Get phi_new using K_inv and current fission source
            phi_new = np.linalg.solve(K, F_phi)

            # Estimate k
            k_new = np.sum(phi_new) / np.sum(phi)
            k_vals.append(k_new)

            # Normalize phi_new
            phi_new /= np.sum(phi_new * h)

            # Compute the L2 norm of the difference
            L2_diff = np.linalg.norm(phi_new - phi, 2)
            L2_diffs.append(L2_diff)

            # Check for convergence
            if len(L2_diffs) > 1 and L2_diffs[-1] < tolerance and abs(
                    1 / k_vals[-1] * (k_vals[-1] - k_vals[-2])) < tolerance:
                break

            # Update phi for next iteration
            phi = phi_new

        # Output
        print('\n*********************** Output for N =', n, '***********************\n')
        print(f'[Number of Iterations]: {len(k_vals)}')
        print(f'[Number of Discrete Spatial bins]: {n}')
        print('\n[Boundary Conditions]:\n(Left) Reflective\n(Right) Zero-Flux')
        print(f'\n[Final k-value after {len(k_vals)} iterations]: {k_vals[-1]:.5f}')
        print(f'[Analytical k-value]: {k_analytical:.5f}')
        print('\n[Numerical k-values for each iteration]:', k_vals)
        print(f'\n[Numerical flux φ [neutrons⋅s⁻¹⋅cm⁻²] at (x = 0)]: {phi[0]:.5f}')
        print(f'[Numerical flux φ [neutrons⋅s⁻¹⋅cm⁻²] at (x=5.5)]: {phi[-1]:.5f}')
        print(f'\n[Final φ {n}-vector]:\n', phi)
        print('\n[The sum of neutrons from x=0 too x=5.5 for all times in the system]:', np.sum(phi * h))
        # Plotting
        plt.plot(x, phi, label='Numerical')
        plt.plot(x, phi_analytical, '--', label='Analytical')
        plt.xlabel('x (cm)')
        plt.ylabel('φ(x) [neutrons⋅s⁻¹⋅cm⁻²] ')
        plt.legend()
        plt.title(f'Flux Eigenfunctions for N:{n}')
        plt.show()

        # Return x, phi, and k for this N
        return x, phi, k_vals[-1]


    print('\n[Problem 2]:')
    # Constants and Input Parameters
    Sigma_t = 0.83  # cm^-1
    Sigma_a = 0.1  # cm^-1
    vSigma_f = 0.125  # cm^-1
    xB = 5.5  # cm
    tolerance = 1e-6
    D = 1 / (3 * Sigma_t)
    h=xB/40

    x = np.linspace(0, xB, 40)  # Adjusted linspace to go from 0 to xB

    # Neutron speed
    v = 1.0  # cm/s

    # Analytical solution
    k_analytical = (4 * xB ** 2 * vSigma_f) / (D * np.pi ** 2 + 4 * Sigma_a * xB ** 2)
    a = np.sqrt((vSigma_f - k_analytical * Sigma_a) / (k_analytical * D))
    phi_analytical = np.cos(a * x)

    h_40 = xB / 40
    h_80 = xB / 80

    # For N=40
    x_40, phi_num_40, k_num_40 = compute_flux_eigenvalue(40)
    ek_40 = abs(k_analytical - k_num_40) / k_analytical
    ephi_40 = np.sqrt(np.sum((phi_analytical - phi_num_40)**2 * h_40) / np.sum(phi_analytical**2 * h_40))

    # For N=80
    x_80, phi_num_80, k_num_80 = compute_flux_eigenvalue(80)
    ek_80 = abs(k_analytical - k_num_80) / k_analytical
    x_80_full = np.linspace(0, xB, 80)
    phi_analytical_80 = np.cos(a * x_80_full)

    # Now, when computing the error for N=80:
    ephi_80 = np.sqrt(np.sum((phi_analytical_80 - phi_num_80)**2 * (h_80) / np.sum(phi_analytical_80**2 * (h_80))))


    order_accuracy_k = np.log(ek_40 / ek_80) / np.log(2)
    order_accuracy_phi = np.log(ephi_40 / ephi_80) / np.log(2)

    print('\n\n*********************** REQUESTED ITEMS for N=40 & N=80 ***********************\n')
    print(f'Error for k(N=40) = {ek_40:.5f}')
    print(f'Error for k(N=80) = {ek_80:.5f}')
    print(f'\nError for φ(N=40) = {ephi_40:.5f} [neutrons⋅s⁻¹⋅cm⁻²]')
    print(f'Error for φ(N=80) = {ephi_80:.5f} [neutrons⋅s⁻¹⋅cm⁻²]')
    print(f'\nOrder accuracy for k: {order_accuracy_k:.5f}')
    print(f'Order accuracy for φ: {order_accuracy_phi:.5f}')
    print('\nEvidently, this method is not well suited for solving this system, as the order-accuracy diminishes with finer resolution (more spatial bins)\n')
    print('\nEverything below is unrequested for problem 2, but was done to see, roughly, the optimal number of spatial bins to minimize the error for k relative to the analytical value.')
    print('\n*********************** (UNREQUESTED) Trying again with N=20 ***********************\n')
    h_20=xB/20
    # For N=20
    x_20, phi_num_20, k_num_20 = compute_flux_eigenvalue(20)
    ek_20 = abs(k_analytical - k_num_20) / k_analytical
    x_20_full = np.linspace(0, xB, 20)
    phi_analytical_20 = np.cos(a * x_20_full)
    # Now, when computing the error for N=20:
    ephi_20 = np.sqrt(np.sum((phi_analytical_20 - phi_num_20)**2 * (0.5) * (h_20) / np.sum(phi_analytical_20**2 * (0.5) * (h_20))))
    print('\n*********************** (UNREQUESTED) Output for N=20 ***********************\n')
    print(f'Error for k(N=20) = {ek_20:.5f}')
    print(f'Error for φ(N=20) = {ephi_20:.5f} [neutrons⋅s⁻¹⋅cm⁻²]')
    print('\n*********************** (UNREQUESTED) Trying again with N=10 ***********************')
    h_10=xB/10
    # For N=10
    x_10, phi_num_10, k_num_10 = compute_flux_eigenvalue(10)
    ek_10 = abs(k_analytical - k_num_10) / k_analytical
    x_10_full = np.linspace(0, xB, 10)
    phi_analytical_10 = np.cos(a * x_10_full)
    # Now, when computing the error for N=2:
    ephi_10 = np.sqrt(np.sum((phi_analytical_10 - phi_num_10)**2 * (0.5) * (h_10) / np.sum(phi_analytical_10**2 * (0.5) * (h_10))))
    print('\n*********************** (UNREQUESTED) Output for N=10 ***********************\n')
    print(f'Error for k(N=10) = {ek_10:.5f}')
    print(f'Error for φ(N=10) = {ephi_10:.5f} [neutrons⋅s⁻¹⋅cm⁻²]')

    print('\nThe accuracy appears to maximize somewhere around 20 groups')

sys.stdout = original_stdout  # Revert back to original stdout

print("Output written to Problem_2_output.txt")

##########################################  PROBLEM 3  ##########################################################################


# Redirecting Output to output file
original_stdout = sys.stdout

with open('Problem_3_output.txt', 'w', encoding='utf-8') as theSack:
    sys.stdout = theSack  # Redirect the output to the file
    # Constants and Input Parameters for Region 1
    Sigma_t1 = 0.83  # cm^-1
    Sigma_a1 = 0.1  # cm^-1
    vSigma_f1 = 0.125  # cm^-1
    D1 = 1 / (3 * Sigma_t1)

    # Constants and Input Parameters for Region 2
    Sigma_t2 = 0.73  # cm^-1
    Sigma_a2 = 0.15  # cm^-1
    vSigma_f2 = 0.115  # cm^-1
    D2 = 1 / (3 * Sigma_t2)

    # Define boundary for the two regions
    x1B = 3  # cm
    x2B = 3.7  # Remaining part for region 2
    xB = x1B + x2B

    # Determine how many points fall in each region
    n = 50
    n1 = round(n * x1B / xB)
    n2 = n - n1

    # Cross sections arrays for the two regions
    Sigma_ts = np.concatenate([Sigma_t1 * np.ones(n1), Sigma_t2 * np.ones(n2)])
    Sigma_as = np.concatenate([Sigma_a1 * np.ones(n1), Sigma_a2 * np.ones(n2)])
    vSigma_fs = np.concatenate([vSigma_f1 * np.ones(n1), vSigma_f2 * np.ones(n2)])
    Ds = np.concatenate([D1 * np.ones(n1), D2 * np.ones(n2)])

    # Neutron speed
    v = 1.0  # cm/s
    xB = 5.5  # cm
    h = xB / n
    tolerance = 10e-8

    x = np.linspace(0, xB, n)

    # Matrix K (Diffusion and Absorption)
    K = np.diag((2 * Ds / h ** 2 + Sigma_as) * np.ones(n))
    K += np.diag((-Ds[:-1] / h ** 2) * np.ones(n - 1), -1)
    K += np.diag((-Ds[1:] / h ** 2) * np.ones(n - 1), 1)
    F = np.diag(vSigma_fs * np.ones(n))

    # Reflective boundary at x = 0
    K[0, 0] = 2 * Ds[0] / h + Sigma_as[0]
    K[0, 1] = -2 * Ds[0] / h

    # Marshak boundary at x = xB
    J_in = 0.0  # Set this to your desired incoming current value
    K[-1, -2] = -Ds[-1] / h
    K[-1, -1] = Ds[-1] / h + 2 * Ds[-1] / h ** 2
    F[-1] = 2 * J_in / h

    # Computing K inverse
    K_inv = np.linalg.inv(K)

    # Analytical a-value for characteristic shape
    k_analytical = (4 * xB ** 2 * vSigma_f1) / (D1 * np.pi ** 2 + 4 * Sigma_a1 * xB ** 2)
    a = np.sqrt((vSigma_f1 - k_analytical * Sigma_a1) / (k_analytical * D1))

    # Power Iteration
    phi = np.ones_like(x)  # Starting guess using a flat profile might be more robust
    phi /= np.sum(phi * h)  # Normalize initial flux to represent a single neutron in the system

    k_vals = []  # To keep track of k-values through iterations
    L2_diffs = []  # To store the L2 norms of the differences

    while True:
        # Compute fission source term
        F_phi = vSigma_fs * phi

        # Get phi_new using K_inv and current fission source
        phi_new = np.linalg.solve(K, F_phi)

        # Estimate k
        k_new = np.sum(phi_new) / np.sum(phi)
        k_vals.append(k_new)

        # Normalize phi_new
        phi_new /= np.sum(phi_new * h)

        # Compute the L2 norm of the difference
        L2_diff = np.linalg.norm(phi_new - phi, 2)
        L2_diffs.append(L2_diff)

        # Check for convergence
        if len(L2_diffs) > 1 and L2_diffs[-1] < tolerance and abs(1 / k_vals[-1] * (k_vals[-1] - k_vals[-2])) < tolerance:
            break

        # Update phi for next iteration
        phi = phi_new

    # Compute the extrapolation distances
    d_left = (2 / 3) * D1 / Sigma_a1
    d_right = (2 / 3) * D2 / Sigma_a2

    # Output
    print(f'\n******************** Problem 3 ********************\n\n[Number of Iterations]: {len(k_vals)}\n[Number of Discrete Spatial Regions]: {n}\n')
    print('[Boundary Conditions]:\n(Left) Reflective\n(Right) Marshak\n')
    print(f'[Final k-value after {len(k_vals)} iterations]: {k_vals[-1]:.5f}')
    print('\n[k-values for each iteration]:', k_vals, '\n')
    print(f'[Numerical flux φ [neutrons⋅s⁻¹⋅cm⁻²] at left boundary (x = 0)]: {phi[0]:.5f}')
    print(f'[Numerical flux φ [neutrons⋅s⁻¹⋅cm⁻²] at right boundary (x=5.5)]: {phi[-1]:.5f}\n')
    print(f'[Final φ {n}-vector]:\n', phi)
    print('\n[The sum of all neutrons from x=0 to x=5.5 for all times]:', np.sum(phi * h))
    # Print the extrapolation distance
    print(f"\n[Extrapolation distance for the right boundary (Region 2)]: {d_right:.5f} cm")

# Completing Output Redirection
sys.stdout = original_stdout  # Reset stdout to console
print('Output written to Problem_3_output.txt')

# Plotting if software is installed
plt.plot(x, phi, label='Numerical')
plt.xlabel('x (cm)')
plt.ylabel('φ(x) [neutrons⋅s⁻¹⋅cm⁻²] ')
plt.legend()
plt.title(f'Flux Eigenfunctions for N: {n}')
plt.show()