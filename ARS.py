from model_BKZ import *
from MSIS_security import *
from MLWE_security import *
import math
from Tools import *

N = 2**8  # ring size

class ParameterSet(object):
    def __init__(self, d, d2, q, kappa, eta, eta2, xi, tau, m, k, n, m2, alpha=1, alpha1=1.2, alpha2=1.2, mode="bimodal"):
        """
        This class represents a parameter set for the signature.
        """
        self.mode = mode

        self.d = d  # ring dimension

        self.d2 = d2  # ring size for N = d * d2^nu
        self.nu = int(math.log(N/d, self.d2))

        self.q = q  # modulus
        self.kappa = kappa  # repetition times for soundness
        self.eta = eta  # infinity norm of s
        self.eta2 = eta2  # infinity norm of r
        self.xi = xi  # infinity norm on c
        self.tau = tau  # constant defined in c
        self.m = m  # Height of A
        self.k = k  # Width of A
        self.n = n  # Height of A_1 and A_2
        self.m1 = self.nu + 2  # length of s_1
        self.m2 = m2  # length of randomness r

        if mode == "bimodal":
            # parameters for rejection sample
            self.alpha = alpha
            self.alpha1 = alpha1
            self.alpha2 = alpha2

            self.s = self.alpha * self.tau * self.eta * math.sqrt(self.k * self.d)  # standard deviation s
            self.B = self.s * math.sqrt(2 * self.m * self.d)  # ||z|| < B
            self.norm = self.B / self.tau + 2 * self.tau  # norm of the vector [s | bc]

            self.Bs1 = math.sqrt(self.nu + 2)  # norm of the vector s_1
            self.s1 = alpha1 * tau * self.Bs1  # standard deviation s_1
            self.s2 = alpha2 * tau * eta2 * math.sqrt(m2 * d)  # standard deviation s_2

            self.B1 = self.s1 * math.sqrt(2 * (self.m1 + 1) * d)  # ||z_1|| < B1
            self.B2 = self.s2 * math.sqrt(2 * m2 * d)  # ||z_2|| < B2
            self.norm2 = 8 * tau * math.sqrt(self.B1 * self.B1 + self.B2 * self.B2)

        elif mode == "convolved":
            S = rot(create_column_vector(d, eta, k))
            max_svalue = max_singular_value(S)
            self.s = math.sqrt(2 * math.log(d - 1 + 2 * d / 0.5) / pi)
            self.sigma = math.sqrt(8) * max_svalue * self.s
            self.B = 1.01 * math.sqrt(d * k) * self.sigma
            self.norm = self.B / tau + 2 * tau

            S1 = rot(create_unit_vector_group(d, self.m1))
            max_svalue1 = max_singular_value(S1)
            self.sigma1 = math.sqrt(8) * max_svalue1 * self.s
            self.B1 = 1.01 * math.sqrt(d * self.m1) * self.sigma1

            S2 = rot(create_column_vector(d, eta2, m2))
            max_svalue2 = max_singular_value(S2)
            self.sigma2 = math.sqrt(8) * max_svalue2 * self.s
            self.B2 = 1.01 * math.sqrt(d * m2) * self.sigma2
            self.norm2 = 8 * tau * math.sqrt(self.B1 * self.B1 + self.B2 * self.B2)

        else:
            print("mode error")


def cal_size(dps):
    print("")
    print("-----------------Size of Signature---------------------")

    log_q = int(math.log(dps.q, 2))

    # full size parameters
    size_c = dps.d * math.ceil(math.log(2 * dps.xi + 1))
    size_full_elems = (dps.n + dps.nu + 2 * dps.kappa) * dps.d * log_q + (2 + dps.k + dps.m * dps.k + dps.m) * dps.d * (log_q + 1)

    # Gaussian size parameters
    size_z = 0
    if dps.mode == "bimodal":
        size_z = dps.k * dps.d * (2.57 + math.ceil(math.log(dps.s, 2))) + dps.m1 * dps.d * (
                    2.57 + math.ceil(math.log(dps.s1, 2))) + dps.m2 * dps.d * (
                         2.57 + math.ceil(math.log(dps.s2, 2)))
    elif dps.mode == "convolved":
        size_z = dps.k * dps.d * (2.57 + math.ceil(math.log(dps.sigma, 2))) + dps.m1 * dps.d * (
                    2.57 + math.ceil(math.log(dps.sigma1, 2))) + dps.m2 * dps.d * (
                         2.57 + math.ceil(math.log(dps.sigma2, 2)))
    else:
        print("mode error")

    # total size of the signature
    total_size = (size_c + size_full_elems + size_z) / 1024 / 8

    print("size:", total_size, "KB, using", dps.mode, "Gaussian")



# MSIS Problem: [A | qj] * [s | bc]^T = 0
def MSIS_1(dps):
    ps = MSISParameterSet(dps.d, dps.k, dps.m, dps.norm, dps.q, "l2")
    x = MSIS_summarize_attacks(ps)
    print("delta: ", delta_BKZ(x[0]))  # hermit factor
    print("svp_classical: ", svp_classical(x[0]))
    print("svp_quantum: ", svp_quantum(x[0]))


# MSIS Problem: [A1 | A2] [cz1-c'z1 | cz2-c'z2]^T = 0
def MSIS_2(dps):
    ps = MSISParameterSet(dps.d, (dps.m1+dps.m2+1), dps.n, dps.norm2, dps.q, "l2")
    x = MSIS_summarize_attacks(ps)
    print("delta: ", delta_BKZ(x[0]))  # hermit factor
    print("svp_classical: ", svp_classical(x[0]))
    print("svp_quantum: ", svp_quantum(x[0]))


# MLWE Problem: b = A_0 s_1 + s_2
def MLWE_1(dps):
    l = dps.k - dps.m - 1
    ps = MLWEParameterSet(dps.d, dps.k, dps.m, dps.eta, dps.q, "uniform")
    x = MLWE_summarize_attacks(ps)
    print("delta: ", delta_BKZ(x[0]))  # hermit factor
    print("svp_classical: ", svp_classical(x[0]))
    print("svp_quantum: ", svp_quantum(x[0]))


# MLWE Problem: [t_A  t_s  t_x  t_g  t]^T = r + [A_1*s_1 s x g f_1]^T
def MLWE_2(dps):
    ps = MLWEParameterSet(dps.d, dps.n+dps.k+dps.m*dps.k+dps.m+dps.nu+dps.k+1, dps.m2, dps.eta2, dps.q, "uniform")
    x = MLWE_summarize_attacks(ps)
    print("delta: ", delta_BKZ(x[0]))  # hermit factor
    print("svp_classical: ", svp_classical(x[0]))
    print("svp_quantum: ", svp_quantum(x[0]))


def security_test(dps):
    print("-----------------Security Test---------------------")
    print("mode:", dps.mode)
    print("")

    print("[MSIS problem in signing]")
    MSIS_1(dps)

    print("")
    print("[MSIS problem in NIZK proof]")
    MSIS_2(dps)

    print("")
    print("[MLWE problem in signing]")
    MLWE_1(dps)

    # print("")
    # print("[MLWE problem in NIZK proof]")
    # MLWE_2(dps)



if __name__ == '__main__':
    """
      using bimodal Gaussian
    """
    dps_96_bimodal = ParameterSet(d=64, d2=4, q=2 ** 29, kappa=10, eta=2, eta2=1, xi=8, tau=140, m=2, k=10, n=17, m2=25,
                                  alpha=1, alpha1=1.2, alpha2=1.2, mode="bimodal")  # parameters for 96 bits security
    # # parameters for 128bits using bimodal Gaussian
    # dps_128_bimodal = ParameterSet(d=64, d2=4, q=2 ** 31, kappa=10, eta=5, eta2=1, xi=8, tau=140, m=3, k=12, n=19,
    #                                m2=25, alpha=1, alpha1=1.2, alpha2=1.2, mode="bimodal")
    #
    # # parameters for 192bits using bimodal Gaussian
    # dps_192_bimodal = ParameterSet(d=128, d2=4, q=2 ** 32, kappa=10, eta=6, eta2=1, xi=2, tau=59, m=2, k=8, n=12,
    #                                m2=25, alpha=1, alpha1=1.2, alpha2=1.2, mode="bimodal")
    security_test(dps_96_bimodal)  # security test
    cal_size(dps_96_bimodal)  # calculate size of signature

    print("")

    """
      using convolved Gaussian
    """
    dps_96_convolved = ParameterSet(d=64, d2=4, q=2 ** 26, kappa=10, eta=8, eta2=1, xi=8, tau=140, m=2, k=7, n=15,
                                    m2=25, mode="convolved")  # parameters for 96 bits security
    security_test(dps_96_convolved)  # security test
    cal_size(dps_96_convolved)  # calculate size of signature








