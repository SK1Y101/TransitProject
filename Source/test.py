import TransitProject

from urllib.request import urlopen
import json

exoclock_planets = json.loads(urlopen('https://www.exoclock.space/database/planets_json').read())
print(exoclock_planets)
