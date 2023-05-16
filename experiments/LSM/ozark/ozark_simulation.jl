t0 = FT(0)
N_days = 60
tf = t0 + FT(3600 * 24 * N_days)
dt = FT(30);
halfhourly = Array(t0:(1800):(t0 + N_days * 3600 * 24))
timestepper = RK4()
