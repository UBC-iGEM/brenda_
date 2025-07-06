# %%
# Imports
from matplotlib import pyplot as plt
from brendapyrser import BRENDA
from numpy import average

# Initialize BRENDA parser
# brenda = BRENDA("brenda_db.txt")

# %%
# Retrieve carbonic anhydrase reaction (EC 4.2.1.1)
rxn = brenda.reactions.get_by_id("4.2.1.1")

# Plot Km values
km_vals = rxn.KMvalues.get_values()
km_vals = [v for v in km_vals if isinstance(v, (int, float))]
km_vals = [v for v in km_vals if 0 <= v <= 1]
print(average(km_vals))

plt.hist(km_vals, bins=30)
plt.title("Km values for EC 4.2.1.1")
plt.xlabel(f"Km ({brenda.units['KM']})")
plt.ylabel("Frequency")
plt.savefig("fig_km.jpg")
plt.clf()

# %%
# Plot pH optima
ph_conditions = rxn.PH
ph_optima_dicts = ph_conditions["optimum"]
ph_optima = [
    entry["value"]
    for entry in ph_optima_dicts
    if isinstance(entry["value"], (int, float))
]
ph_optima = [v for v in ph_optima if 0 <= v <= 14]
print(average(ph_optima))

plt.hist(ph_optima, bins=14)
plt.title("pH optima for EC 4.2.1.1")
plt.xlabel("pH")
plt.ylabel("Frequency")
plt.savefig("fig_ph.jpg")
plt.clf()

# %%
# Plot pH ranges
ph_range_dicts = ph_conditions["range"]
ph_ranges = [
    entry["value"]
    for entry in ph_range_dicts
    if (
        isinstance(entry["value"], list)
        and len(entry["value"]) == 2
        and all(isinstance(v, (int, float)) for v in entry["value"])
    )
]
fig, ax = plt.subplots()
for i, (low, high) in enumerate(ph_ranges):
    ax.hlines(y=i, xmin=low, xmax=high, linewidth=2)
ax.set_xlabel("pH")
ax.set_ylabel("Sample index")
ax.set_title("pH measurement ranges for EC 4.2.1.1")
ax.set_xlim(0, 14)
plt.savefig("fig_ph_ranges.jpg")
plt.clf()

# %%
# Plot temperature optima
temp_conditions = rxn.temperature
temp_optima_dicts = temp_conditions["optimum"]
temp_optima = [
    entry["value"]
    for entry in temp_optima_dicts
    if isinstance(entry["value"], (int, float))
    if -200 <= entry["value"] <= 200
]
print(average(temp_optima))

plt.hist(temp_optima, bins=30)
plt.title("Temperature optima for EC 4.2.1.1")
plt.xlabel("Temp (deg C)")
plt.ylabel("Frequency")
plt.savefig("fig_temp.jpg")
plt.clf()

# %%
# Plot temp ranges
temp_range_dicts = temp_conditions["range"]
temp_ranges = [
    entry["value"]
    for entry in temp_range_dicts
    if (
        isinstance(entry["value"], list)
        and len(entry["value"]) == 2
        and all(isinstance(v, (int, float)) for v in entry["value"])
        and all(-200 <= v <= 200 for v in entry["value"])
    )
]
fig, ax = plt.subplots()
for i, (low, high) in enumerate(temp_ranges):
    ax.hlines(y=i, xmin=low, xmax=high, linewidth=2)
ax.set_xlabel("Temp (deg C)")
ax.set_ylabel("Sample index")
ax.set_title("Temp ranges for EC 4.2.1.1")
plt.savefig("fig_temp_ranges.jpg")
plt.clf()

# %%
# List data
data = {
    "organism": rxn.organisms,
    "ph_optima": [
        (val if isinstance(val, (int, float)) and 0 <= val <= 14 else None)
        for entry in rxn.PH.get("optimum", [])
        for val in (entry.get("value"),)
    ],
    "temp_optima": [
        (val if isinstance(val, (int, float)) and -200 <= val <= 200 else None)
        for entry in rxn.temperature.get("optimum", [])
        for val in (entry.get("value"),)
    ],
}

print(len(data["organism"]))
print(len(data["ph_optima"]))
print(len(data["temp_optima"]))

records = [
    {"organism": org, "ph_optima": ph, "temp_optima": temp}
    for org, ph, temp in zip(
        data["organism"], data["ph_optima"], data["temp_optima"]
    )
]

print(records)
