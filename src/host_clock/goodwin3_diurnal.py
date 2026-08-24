import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# real 6-point profile (ZT, mmol/gDW/h) -- from your butyrate_timecourse.py run
zt = np.array([1, 5, 9, 13, 17, 21])
vals = np.array([2.9331, 2.8431, 2.7699, 2.6228, 3.1592, 3.5366])

zt_ext = np.concatenate([zt, [zt[0] + 24]])
vals_ext = np.concatenate([vals, [vals[0]]])
B_spline = CubicSpline(zt_ext, vals_ext, bc_type="periodic")

def B_of_t(t):
    return float(B_spline(t % 24))

VS, KS1, KS2, KD, KI, N = 1.0, 1.0, 1.0, 0.154, 1.0, 8
ALPHA = 0.15

def rhs_diurnal(t, y):
    X, Y, Z = y
    B = B_of_t(t)
    dX = VS * (KI**N / (KI**N + Z**N)) * (1 + ALPHA * B) - KD * X
    dY = KS1 * X - KD * Y
    dZ = KS2 * Y - KD * Z
    return [dX, dY, dZ]

def rhs_baseline(t, y):
    X, Y, Z = y
    dX = VS * (KI**N / (KI**N + Z**N)) * (1 + ALPHA * 0.0) - KD * X
    dY = KS1 * X - KD * Y
    dZ = KS2 * Y - KD * Z
    return [dX, dY, dZ]

T = 480
t_eval = np.linspace(0, T, T*20)
sol_base = solve_ivp(rhs_baseline, (0,T), [0.3,0.3,0.3], t_eval=t_eval, method="RK45", rtol=1e-9, atol=1e-9)
sol_diur = solve_ivp(rhs_diurnal, (0,T), [0.3,0.3,0.3], t_eval=t_eval, method="RK45", rtol=1e-9, atol=1e-9)

def analyze(sol, label):
    mask = sol.t >= T-120
    x = sol.y[0][mask]; t = sol.t[mask]
    peaks = [t[i] for i in range(1,len(x)-1) if x[i]>x[i-1] and x[i]>x[i+1]]
    period = np.mean(np.diff(peaks)) if len(peaks)>2 else float("nan")
    amp = (x.max()-x.min())/2
    mean = x.mean()
    print(f"{label}: period={period:.3f}h, amplitude={amp:.5f}, mean={mean:.5f}, n_peaks_48h={len(peaks)}")
    return period, amp, mean, peaks

print("="*65)
print("Host clock coupled to REAL 6-point diurnal butyrate profile")
print("="*65)
pb = analyze(sol_base, "Baseline (B=0)")
pd_ = analyze(sol_diur, "Diurnal B(t), real 6-point profile")

if len(pb[3])>1 and len(pd_[3])>1:
    base_phase = pb[3][-1] % 24
    diur_phase = pd_[3][-1] % 24
    print(f"\nBaseline last peak time (mod 24h): {base_phase:.2f}")
    print(f"Diurnal last peak time (mod 24h): {diur_phase:.2f}")
    print(f"Amplitude change (diurnal vs baseline): {(pd_[1]-pb[1])/pb[1]*100:.1f}%")
    print(f"Mean expression change: {(pd_[2]-pb[2])/pb[2]*100:.1f}%")

print(f"\nPeriod under real diurnal forcing: {pd_[0]:.3f}h (baseline free-running: {pb[0]:.3f}h)")

fig, axes = plt.subplots(2, 1, figsize=(9.5, 7), sharex=True)
mask_plot = sol_base.t >= T-96
tp = sol_base.t[mask_plot] - (T-96)

axes[0].plot(tp, sol_base.y[0][mask_plot], label="Per mRNA (baseline, B=0)", color="#4C72B0")
axes[0].plot(tp, sol_diur.y[0][mask_plot], label="Per mRNA (real diurnal butyrate B(t))", color="#C44E52")
axes[0].set_ylabel("Per mRNA (a.u.)")
axes[0].legend()
axes[0].set_title("Host clock coupled to REAL 6-point diurnal butyrate export profile")

t_full = np.linspace(0, 96, 2000)
axes[1].plot(t_full, [B_of_t(tt) for tt in t_full], color="#55A868", label="B(t) (periodic spline through 6 real points)")
for i in range(4):
    axes[1].scatter(zt + i*24, vals, color="#333333", zorder=5, s=25,
                     label="Real optimized ZT points" if i==0 else None)
axes[1].set_ylabel("Butyrate export capacity\n(mmol/gDW/h)")
axes[1].set_xlabel("Time (h)")
axes[1].legend()

plt.tight_layout()
plt.savefig("/home/aceren/diurnal_host_microbiome/diurnal_coupling.png", dpi=150)
print("\nSaved: /home/aceren/diurnal_host_microbiome/diurnal_coupling.png")
