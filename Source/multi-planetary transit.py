# import modules
from astropy import units as u

# import my module
from ProjectModules import nbody

# initialise our system
def main():
    sun = nbody.Body(1.98855*u.kg)

    earth = nbody.Body(5.72E24 * u.kg)
    earth.orbit(sun, 1 * u.au, 0, 0)

# execute if we are the top level code
if __name__ == "__main__":
    main()
