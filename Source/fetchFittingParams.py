import pandas as pd
import numpy as np

df = pd.read_csv("TTVTestModelFittingParameters.csv")

p1 = pd.DataFrame(columns=["Model", "Method", "mu", "a", "e", "omega"])
p2 = pd.DataFrame(columns=["Model", "Method", "mu", "a", "e", "omega"])
p3 = pd.DataFrame(columns=["Model", "Method", "mu", "a", "e", "omega"])

for idx in df.index:
    thisfit = df.iloc[idx]

    name = thisfit["models"]
    params = [x.strip() for x in thisfit["solutions"].replace("[", "").replace("]", "").split(" ") if x][:-1]

    per = int(name.split(":")[1].split(",")[0].replace("perturbers", "").strip())

    model = name.split(":")[0].replace("model", "Model ")
    method = name.split(",")[1].replace("differential", "diff.").replace("evolution", "evo.").replace("annealing", "an.").replace("fitting", "")

    for i, x in enumerate(np.reshape(params, (-1, 6))):
        mu = round(float(x[1])*10**6, 3)
        a = round(float(x[0])*10**6, 3)
        e = round(float(x[2]), 3)
        omega = round(float(x[4]), 3)

        thisrow = {"Model":"" if i else model, "Method":"" if i else method, "mu":mu, "a":a, "e":e, "omega":omega}

        if per == 1:
            p1 = p1.append(thisrow, ignore_index=True)
        if per == 2:
            p2 = p2.append(thisrow, ignore_index=True)
        if per == 3:
            p3 = p3.append(thisrow, ignore_index=True)

print(p1.to_latex(index=False))
print(p2.to_latex(index=False))
print(p3.to_latex(index=False))
