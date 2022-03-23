from decimal import *

D = Decimal

Fixed = Decimal(10) ** -30

def tofixed(x):
    return D(x).quantize(Fixed)

# taken from deminal library
def pi():
    pi_ = D(3)
    _pi = 0
    t = D(3)
    n = 1
    na = 0
    d = 0
    da = 24
    while _pi != pi_:
        _pi = pi_
        n, na = n+na, na+8
        d, da = d+da, da+32
        t = (t*n)/d
        pi_ += t
    return +pi_

def sin(x):
    # Accurate sine using taylor expansion
    x = D(x)
    #increase precision
    # sin(x) ~ x - x**3 / 3! + x**5 / 5! - x**7 / 7! + ...
    _sinx, sinx = 0, x
    p, f, s = 1, 1, 1
    xp = x
    # keep going until we match up to our precision
    while sinx != _sinx:
        _sinx = sinx
        # compute the new power (sin is zero for evey other power)
        p += 2
        # compute the new factorial
        f *= p * (p-1)
        # compute our x to the power
        xp *= x * x
        # compute the new sign
        s *= -1
        # add this term to our sine value
        sinx += s * xp / f
    # decrease before returning
    return +sinx

def arcsin(x):
    x = D(x)
    # the taylor series is slow to converge near x=1, so we'll reduce the range
    if x>D(0.7071):
        return D(0.5)*pi() - arcsin((1 - x**2).sqrt())
    asinx, asinx_ = x, 0
    #asin(x) ~ x + (1/2)x^3 / 3 + ...
    p = 1
    f = 1
    xp = x
    i=0
    while asinx != asinx_:
        asinx_ = asinx
        # power
        p += 2
        # fraction
        f *= D(p-2) / D(p-1)
        # power of x
        xp *= x * x
        # new approximation
        asinx += xp * f / D(p)
    return +asinx

def cos(x):
    x = D(x)
    # cos(x) ~ 1 - x**2 / 2! + ...
    _cosx, cosx = 0, 1
    p, f, s = 0, 1, 1
    xp = 1
    # keep going until we match up to our precision
    while cosx != _cosx:
        _cosx = cosx
        p +=2
        xp *= x*x
        f *= p * (p-1)
        s*=-1
        cosx += s * xp / f
    # decrease before returning
    return +cosx

def arccos(x):
    x = D(x)
    pi_ = pi()
    asinx = arcsin(x)
    acosx = pi_ * D(0.5) - asinx
    return +acosx

def tan(x):
    x = D(x)
    tanx = sin(x) / cos(x)
    return +tanx

def arctan(x):
    x = D(x)
    # arctan(x) ~ x - x**3 / 3 + x**5 / 5 - ...
    atanx, _atanx = x, 0
    xp = x
    p = 1
    s = 1
    while atanx != _atanx:
        _atanx = atanx
        p += 2
        xp *= x * x
        s *= -1
        atanx += s * xp / p
    # decrease before returning
    return +atanx

def arctan2(y, x):
    x, y = D(x), D(y)
    # if x is zero
    if x == 0 and y>0:
        atan = pi() * D(0.5)
    elif x == 0 and y < 0:
        atan = -pi() * D(0.5)
    # otherwise, compute atan
    else:
        atan = arctan(y / x)
    # compute additional stuff
    if x < 0 and y >= 0:
        atan += pi()
    elif x < 0 and y < 0:
        atan -= pi()
    return +atan
