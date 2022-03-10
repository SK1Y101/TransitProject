

# compute the area of intersection of two circles, given the distance from their centres and their radii
def area(dist, rad1, rad2):
    # convert to Decimal if not already
    dist, rad1, rad2 = D(dist), D(rad1), D(rad2)
    # if they are too far apart, area covered is zero
    if dist > rad1+rad2:
        return D(0)
    # if one is completely covered by the other
    if dist < max(rad1, rad2) - min(rad1, rad2):
        return D(0.5) * trig_funcs.pi() * min(rad1, rad2)**2
    return D(0.5) * trig_funcs.pi() * min(rad1, rad2)**2
    # distance from line between centres to edge of overlap
    dr = D(D(1) / ((D(1) / rad1**2) + (D(1) / rad2**2))).sqrt()
    # distance from the centre of each circle to the centre of the intersection
    d1 = D(rad1**2 - dr**2).sqrt()
    d2 = D(rad2**2 - dr**2).sqrt()
    # angle between the line between the centres and the edge of the overlap
    t1, t2 = trig_funcs.arccos(d1 / rad1), trig_funcs.arccos(d2 / rad2)
    # area of the circle secgement
    a1, a2 = D(0.5) * rad1**2 * t1, D(0.5) * rad2**2 * t2
    # area bounded by the triangular section
    at1, at2 = D(0.5) * d1 * dr, D(0.5) * d2 * dr
    # total area of the intersection is the circle segment - triangle for each circle
    return a1 + a2 - at1 - at2

system.depths = []
system.times = []

lim = (primary_sma * (resonance**(3/2) if resonance > 1 else 1)).value
for body in system.bodies:
    body.size = float(D(225) * (body.radius / D(lim)))**2

def compute_occlusion(star, body):
# compute the angular seperation between the target and star
# vector to observer
p_vec, p_mag = body.dist, body.dist_mag
s_vec, s_mag = star.dist, star.dist_mag

# dot and magnitude product
dot_prod = sum(p_vec * s_vec)
mag_prod = p_mag * s_mag

# angle between
cosphi = dot_prod / mag_prod
phi = D(trig_funcs.arccos(cosphi))

# stop if the planet is not near to occluding the star
if phi > body.theta + star.theta:
    return D(0)

# convert angular values above to distances at the projected point
dist_ang_rad = trig_funcs.tan(phi) * D(star.dist_mag)
#print(dist_ang_rad)
body_ang_rad = trig_funcs.tan(body.theta) * D(star.dist_mag)
#print(body_ang_rad)
star_ang_rad = trig_funcs.tan(star.theta) * D(star.dist_mag)
#print(star_ang_rad)
# compute the area of the stars disk covered by the planet
a = area(dist_ang_rad, body_ang_rad, star_ang_rad)
return a

def simulate_step(keepBarycentreConstant=True):
# propogate one timestep
system.propogate(system.tstep, keepBarycentreConstant=True)

# loop over all planets
for body in system.bodies:
    # compute distance to body
    body.dist = (observer_pos - body.pos)
    body.dist_mag = np.linalg.norm(body.dist)#np.sqrt(sum(body.dist * body.dist))

    # compute the angular radius
    body.theta = trig_funcs.arctan2(body.radius, body.dist_mag)

    if body == star:
        continue
    # if the planet is on the front half of the star
    if body.dist_mag < star.dist_mag:
        # fetch the area of the star covered by this planet
        body.area_covered = compute_occlusion(star, body)
    else:
        body.area_covered = D(0)

# fetch the occlusion areas
area = [body.area_covered for body in system.bodies if body != star]
# if no planet is transiting
if area.count(D(0)) == len(area):
    depth = 1
# if only a single planet is transiting
elif area.count(D(0)) == len(area)-1:
    # fetch the area
    area = [a for a in area if a!= D(0)][0]
    # and angular radius of the star
    star_ang_rad = trig_funcs.tan(star.theta) * D(star.dist_mag)
    # compute the transit depth
    depth = D(1) - min(D(1), area / D(D(0.5) * trig_funcs.pi() * star_ang_rad**2))
else:
    pass
    #print(compute_occlusion(p1, p2))
# combine the transit depths

# compute the difference in transit depth due to the star wobbling
obs_dist = np.linalg.norm(system.barycentre - observer_pos)
extra_dist = np.linalg.norm(system.barycentre - star.pos)
wobble = D( ((obs_dist + extra_dist) / obs_dist) )
depth *= wobble

system.times.append(system.time)
system.depths.append(depth)

# simulate a timestep and plot it
def plotsimulatedtimestep(ax, i, keepBarycentreConstant=True):
simulate_step(keepBarycentreConstant)
#ax[0].scatter(*[float(val) for val in observer_pos], c="cyan")
for body in system.bodies:
    # plot the previous position of the body
    ax[0].plot(*[[float(v) for v in val] for val in system.positions[body].T], c=body.colour, lw=0.5)
    # plot the current position of the body
    ax[0].scatter(*[float(val) for val in body.pos], c=body.colour, s=body.size)
# and plot on the lower graph
ax[1].plot(system.times, system.depths)
ax[0].set_xlim([-lim, lim])
ax[0].set_ylim([-lim, lim])
ax[0].set_zlim([-lim, lim])
#system.plot(updatefunc=plotsimulatedtimestep, interval=10, keepBarycentreConstant=True)
