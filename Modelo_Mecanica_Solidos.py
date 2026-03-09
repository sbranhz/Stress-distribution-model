import numpy as np
import matplotlib.pyplot as plt

# ─── Parámetros de la sección transversal ───────────────────────────────────
ac = 50        # Ancho en la parte superior (mm)
ap = 70        # Ancho en la parte inferior (mm)
h  = 100       # Altura de la viga (mm)
b  = 1         # Espesor (mm)

# ─── Parámetros de carga y material ─────────────────────────────────────────
g   = 9.8e3    # Gravedad (mm/s²)
f   = 1.04e6   # Carga distribuida aplicada (N/mm)
gc  = 2.3 * g  # Peso específico del material (N/mm³)
r   = 10       # Resolución de la malla (puntos)

# ─── Flags de configuración ──────────────────────────────────────────────────
body = True    # Incluir peso propio de la viga
uni  = False   # True = carga uniforme, False = carga triangular

# ─── Peso propio y geometría ─────────────────────────────────────────────────
w = (ac + ap) * h * b * gc / 2
y = np.linspace(h, 0, r)
Y = np.zeros((r, r))
for i in range(r):
    for j in range(r):
        Y[i, j] = y[i]

m = (ac - ap) / h
a = m * y + ap

# ─── Cortante basal ──────────────────────────────────────────────────────────
if uni:
    vb = f * h
else:
    vb = f * h / 2

# ─── Centro de gravedad del peso propio ──────────────────────────────────────
if body:
    cb = (ap / 2) - (ac * (3 * ac + (ap - ac) ** 2) / (3 * ac + 3 * ap))
else:
    cb = 0

# ─── Momento en la base ──────────────────────────────────────────────────────
if uni:
    Mb = vb * h / 2
else:
    Mb = vb * h / 3
if body:
    Mb = Mb - w * cb

# ─── Momento flector a lo largo de y ─────────────────────────────────────────
if uni:
    My = -Mb + f * y ** 2 / 2
else:
    t  = f - f * y / h
    My = -Mb + vb * y - (t + f) * y ** 2 / 4

# ─── Esfuerzo normal Syy ─────────────────────────────────────────────────────
Eny = a / 2
cx  = np.zeros((r, r))
for i in range(r):
    for j in range(r):
        cx[i, j] = a[i] * j / (r - 1)

Iy   = b * a ** 3 / 12
Syy  = np.zeros((r, r))
for i in range(r):
    for j in range(r):
        if body:
            Syy[i, j] = -My[i] * (Eny[i] - cx[i, j]) / Iy[i] - (a[i] + ac) * (h - y[i]) / 2
        else:
            Syy[i, j] = -My[i] * (Eny[i] - cx[i, j]) / Iy[i]

Sy  = np.zeros((r, r))
Sy1 = np.zeros((r, r))
for i in range(r):
    for j in range(int(r / 2)):
        Sy[i,  j + int(r / 2)] = Syy[i, j + int(r / 2)]
        Sy1[i, j]              = Syy[i, j]
Sy1  = np.flip(Sy1, 0)
Syyy = Sy1 + Sy
Syy1 = np.flip(Syy, 1)

# ─── Esfuerzo cortante Txy ───────────────────────────────────────────────────
if uni:
    vy = vb - f * y
else:
    vy = vb - (t + f) * y / 2

Q = np.zeros((r, r))
if r % 2 != 0:
    r2 = int((r + 1) / 2)
    for i in range(r):
        for j in range(r2):
            Q[i, j] = cx[i, j] * (Eny[i] - cx[i, j] + cx[i, j] / 2)
    for i in range(r):
        for j in range(r2 - 1):
            Q[i, j + r2] = Q[i, r2 - 2 - j]
else:
    r2 = int(r / 2)
    for i in range(r):
        for j in range(r2):
            Q[i, j] = cx[i, j] * (Eny[i] - cx[i, j] + cx[i, j] / 2)
    for i in range(r):
        for j in range(r2):
            Q[i, j + r2] = Q[i, r2 - 1 - j]

Txy = np.zeros((r, r))
for i in range(r):
    for j in range(r):
        Txy[i][j] = vy[i] * Q[i, j] / (Iy[i] * b)

# ─── Esfuerzo normal Sxx ─────────────────────────────────────────────────────
x  = np.linspace(ap, 0, r)
X  = np.zeros((r, r))
for i in range(r):
    for j in range(r):
        X[i, j] = x[i]

if ac == ap:
    hx = np.ones(r) * h
else:
    hx    = (x - ap) / m
    hx[0] = 1e-100

l = 6 * Mb / ap ** 2
d = l - 2 * l / ap * x

if body:
    fb = np.zeros(r)
    for i in range(r):
        if x[i] <= ac:
            fb[i] = ap * h * b * gc / (r - 1)
        else:
            fb[i] = (hx[i - 1] + 2 * hx[i] + hx[i + 1]) * ap * b * gc / (4 * (r - 1))
    if ac == ap:
        fb[0] = fb[1]
    fb[0] = hx[1] * ap / (8 * (r - 1))
    d = d - fb

vx    = (l + d) * x / 2
vx[0] = 0

c = np.zeros(r)
for i in range(r):
    if x[i] <= ap / 2:
        c[i] = x[i] * (d[i] + 2 * l) / (3 * (d[i] + l))
    else:
        c[i] = (3 * l * ap * (x[i] - ap / 6) + 2 * abs(d[i]) * (x[i] - ap / 2) ** 2) / \
               (3 * l * ap + 6 * abs(d[i]) * (x[i] - ap / 2))

Enx = np.zeros(r)
for i in range(r):
    if x[i] <= ac:
        Enx[i] = h / 2
    else:
        Enx[i] = hx[i] / 2

Mx = vx * c
if not uni:
    for i in range(r):
        Mx[i] = Mx[i] - (f * h / 2) * (Enx[i] - h / 3)

cy = np.zeros((r, r))
for i in range(r):
    for j in range(r):
        if x[i] <= ac:
            cy[i, j] = h * j / (r - 1)
        else:
            cy[i, j] = hx[i] * j / (r - 1)

Ix = np.zeros(r)
for i in range(r):
    if x[i] <= ac:
        Ix[i] = b * h ** 3 / 12
    else:
        Ix[i] = b * hx[i] ** 3 / 12
Ix[0] = 1

f1 = np.ones(r) * f
Sw = np.zeros(r)
for i in range(r):
    if ac == ap:
        Sw[i] = f1[i]
    else:
        if hx[i] >= h:
            Sw[i] = f1[i]
        else:
            Sw[i] = f * b * h / (b * hx[i])
Sw[0] = f * 1e5

Sxx = np.zeros((r, r))
for i in range(r):
    for j in range(r):
        Sxx[i, j] = Mx[i] * (Enx[i] - cy[i, j]) / Ix[i] - Sw[i]
        if ac != ap:
            Sxx[0, j] = Sxx[1, j] * 2
Sxx = np.flip(np.flip(Sxx, 1), 0)

# ─── Esfuerzos principales ───────────────────────────────────────────────────
Sc   = (Sxx + Syy) / 2
Tmax = (((Sxx - Syy) / 2) ** 2 + Txy ** 2) ** 0.5
Smin = np.flip(Tmax - Sc, 1)
Smax = Tmax + Sc

# ─── Visualización ───────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
fig.suptitle("Distribución de Esfuerzos — Viga Trapezoidal", fontsize=14, fontweight='bold')

plots = [
    (Syy,  cx, Y, "Esfuerzo Normal σyy"),
    (Txy,  cx, Y, "Esfuerzo Cortante τxy"),
    (Sxx,  cx, Y, "Esfuerzo Normal σxx"),
    (Smin, cx, Y, "Esfuerzo Principal Mínimo"),
    (Smax, cx, Y, "Esfuerzo Principal Máximo"),
    (Tmax, cx, Y, "Esfuerzo Cortante Máximo"),
]

for ax, (data, xg, yg, title) in zip(axes.flat, plots):
    cf = ax.contourf(xg, yg, data, 100, cmap='RdYlBu_r')
    fig.colorbar(cf, ax=ax)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")

plt.tight_layout()
plt.savefig("stress_distribution.png", dpi=150, bbox_inches='tight')
plt.show()
print("Análisis completado. Figura guardada como stress_distribution.png")
