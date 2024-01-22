def custom_imf(mass):
    if mass < 0.08:
        return 0
    elif mass < 0.5:
        return 2*mass**(-1.3)
    elif mass < 1:
        return mass**(-2.3)
    elif mass < 150 or mass == 150:
        return mass**(-2.3)
    else:
        return 0
