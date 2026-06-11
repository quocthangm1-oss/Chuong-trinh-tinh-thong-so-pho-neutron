import math


def Q0_alpha(Q0, Er, alpha):
    term1 = (Q0 - 0.429) / (Er ** alpha)
    term2 = 0.429 / ((2 * alpha + 1) * (0.55 ** alpha))
    return term1 + term2


def F_alpha(alpha, a, b, params):

    Q95 = Q0_alpha(params['Q0_95'], params['Er_95'], alpha)
    Q97 = Q0_alpha(params['Q0_97'], params['Er_97'], alpha)
    QAu = Q0_alpha(params['Q0_Au'], params['Er_Au'], alpha)

    termAu = (a - b) * QAu * params['GepAu_GthAu']
    term95 = -a * Q95 * params['Gep95_Gth95']
    term97 = b * Q97 * params['Gep97_Gth97']

    return term95 + term97 + termAu


def solve_alpha_iterative(
        a,
        b,
        params,
        alpha_min=-0.99,
        alpha_max=0.99,
        step=0.00001,
        tol=0.0000000000001):

    alpha = alpha_min
    iterations = 0

    while alpha <= alpha_max:

        Fval = F_alpha(alpha, a, b, params)

        if abs(Fval) <= tol:
            return alpha, iterations + 1, True

        alpha += step
        iterations += 1

    return None, iterations, False


def calc_f(alpha, params):

    Q95 = Q0_alpha(params['Q0_95'], params['Er_95'], alpha)
    Q97 = Q0_alpha(params['Q0_97'], params['Er_97'], alpha)

    numerator = (
        params['Gep95']
        * Q95
        * (params['k0_95'] * params['eps95'])
        / (params['k0_97'] * params['eps97'])
        -
        params['Gep97']
        * Q97
        * (params['Asp95'] / params['Asp97'])
    )

    denominator = (
        params['Gth97']
        * (params['Asp95'] / params['Asp97'])
        -
        params['Gth95']
        * (params['k0_95'] * params['eps95'])
        / (params['k0_97'] * params['eps97'])
    )

    return numerator / denominator


params = {

    # Q0
    'Q0_95': 5.31,
    'Q0_97': 251.6,
    'Q0_Au': 15.7,

    # Er
    'Er_95': 6260,
    'Er_97': 338,
    'Er_Au': 5.65,

    # G-factor
    'Gep95_Gth95': 1.0,
    'Gep97_Gth97': 1.0,
    'GepAu_GthAu': 1.0,

    'Gep95': 1.0,
    'Gep97': 1.0,

    'Gth95': 1.0,
    'Gth97': 1.0,

    # k0
    'k0_95': 0.0000890,
    'k0_97': 0.0000124,

    # hiệu suất
    'eps95': 0.001529,
    'eps97': 0.001488,

    # hoạt độ riêng
    'Asp95': 1.07E5,
    'Asp97': 7.07E4,
}


a = -0.00083809
b = -0.000144899


alpha, n_iter, converged = solve_alpha_iterative(
    a,
    b,
    params,
    step=0.00001,
    tol=0.0000000000001
)

if not converged:
    print("Không hội tụ alpha.")
    quit()


f_factor = calc_f(alpha, params)


Q0_Lu_ref = 1.908933
Er_Lu = 0.158

s0_Au_ref = 17.24
s0_Lu_ref = 1.67

W_Au = 0.055

FCd = 0.991

Gth = 1.0
Gr = 1.0
gTn = 1.0


k0_Au = 1.0
eps_Au = 0.002702194
Asp_Au = 7910611644

k0_Lu = 0.0714
eps_Lu = 0.004854436
Asp_Lu = 802130558


Q0_Au_alpha = Q0_alpha(
    params['Q0_Au'],
    params['Er_Au'],
    alpha
)

Q0_Lu_alpha = Q0_alpha(
    Q0_Lu_ref,
    Er_Lu,
    alpha
)


s0_Au_alpha = s0_Au_ref * (params['Er_Au'] ** (-alpha))
s0_Lu_alpha = s0_Lu_ref * (Er_Lu ** (-alpha))


W_prime_Au = W_Au * (params['Er_Au'] ** (-alpha))


RCd = ((f_factor / Q0_Au_alpha) + 1.0) / FCd


denominator = (
    (RCd * FCd)
    * (
        (gTn / 2.2931)
        / ((1 + 2 * alpha) * (0.55 ** alpha))
        -
        (2 / math.sqrt(math.pi))
        * W_prime_Au
        +
        s0_Au_alpha * Gr
    )
    -
    s0_Au_alpha * Gr
)

R_alpha = (Gth * gTn) / denominator

g_Lu_Tn = (
    (
        (Asp_Lu / (Gth * k0_Lu * eps_Lu))
        /
        (Asp_Au / (k0_Au * eps_Au))
    )
    *
    (Gth * gTn + Gr * R_alpha * s0_Au_alpha)
    -
    Gr * R_alpha * s0_Lu_alpha
) / Gth


print("\n==============================")
print("KẾT QUẢ TÍNH TOÁN")
print("==============================")

print(f"Alpha = {alpha:.6f}")
print(f"Số vòng lặp = {n_iter}")
print(f"f = {f_factor:.6f}")

print("\n----- Q0(alpha) -----")
print(f"Q0_Au(alpha) = {Q0_Au_alpha:.6f}")
print(f"Q0_Lu(alpha) = {Q0_Lu_alpha:.6f}")

print("\n----- S0(alpha) -----")
print(f"S0_Au(alpha) = {s0_Au_alpha:.6f}")
print(f"S0_Lu(alpha) = {s0_Lu_alpha:.6f}")

print("\n----- W'(alpha) -----")
print(f"W'_Au(alpha) = {W_prime_Au:.6f}")

print("\n----- RCd -----")
print(f"RCd = {RCd:.6f}")

print("\n----- r(alpha) -----")
print(f"r(alpha) = {R_alpha:.6f}")

print("\n----- gLu(Tn) -----")
print(f"g_Lu(Tn) = {g_Lu_Tn:.6f}")
