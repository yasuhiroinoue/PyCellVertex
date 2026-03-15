from __future__ import annotations

# Parity constants from C++ `_parameters.h` / `_parameters.cpp`
DELTA_TIME = 1e-4

# reconnection-related params in C++ header
L_THRESHOLD = 1e-1
L_RECONNECTED = 5e-2
TIME_RECONNECT = 10.0
STEP_RECONNECT_DEFAULT = int(TIME_RECONNECT / DELTA_TIME)
