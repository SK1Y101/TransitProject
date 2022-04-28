''' List of TTV models '''

# 2 planet system, zero eccentricity, coplanar, non-resonant.

TTV1 = -sma_1 * mu_1 * p_1 * np.sin(2*np.pi*(n*p_2 - t) / p_1) / (2*np.pi*sma_2)
# where n is the transit number, and t is the time of the transit, such that n*p2 gives the position of the outer planet
