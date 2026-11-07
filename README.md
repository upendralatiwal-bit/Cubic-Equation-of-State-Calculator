# Cubic Equation of State Calculator

A browser-based tool for predicting the PVT behaviour and thermodynamic properties of pure real gases using cubic equations of state.

**Group members:** Upendra Kumar (2025CH11982), Sadiya (2025CH12131)  
Department of Chemical Engineering, IIT Delhi

**Live site:** `https://<your-username>.github.io/cubic-eos/`

## Features
- **EOS models:** van der Waals, Redlich–Kwong, Soave–Redlich–Kwong, Peng–Robinson (generic σ–ε–Ω–Ψ form)
- **PVT:** compressibility factor Z and molar volume V from analytic (Cardano) roots; liquid / vapour / unstable roots with stable-phase selection by minimum Gibbs energy
- **Thermodynamic properties:** residual enthalpy Hᴿ and entropy Sᴿ, plus total H and S relative to an ideal-gas reference state (T₀, P₀) using ideal-gas Cp polynomials
- **Fugacity:** pure-component fugacity coefficient φ and fugacity f = φP
- **Plots:** P–V isotherms (sub-/near-/super-critical), Z vs P, Hᴿ & Sᴿ vs P, fugacity vs P, and Z vs P comparison of all four EOS
- **Database:** 12 common gases (methane, ethane, propane, n-butane, ethylene, CO₂, CO, N₂, O₂, H₂, NH₃, H₂O); Tc, Pc, ω are editable for any other gas

## Key equations
P = RT/(V−b) − a(T)/[(V+εb)(V+σb)]

Hᴿ/RT = Z − 1 + (d ln α/d ln Tr − 1) q I  
Sᴿ/R = ln(Z − B) + (d ln α/d ln Tr) q I  
ln φ = Z − 1 − ln(Z − B) − q I

where B = bP/RT, q = a/(bRT), I = ln[(Z+σB)/(Z+εB)]/(σ−ε) (I = B/Z for vdW).

## Files
| File | Purpose |
|---|---|
| `index.html` | User interface and plots (Chart.js) |
| `eos.js` | EOS engine: root solver, residual properties, fugacity, component data |

## Run locally
Just open `index.html` in a browser. The engine also runs in Node:
```js
const {solve, COMPONENTS} = require("./eos.js");
console.log(solve("PR", COMPONENTS["Methane"], 300, 50).stable);
```

## Reference
Smith, J.M., Van Ness, H.C., Abbott, M.M., *Introduction to Chemical Engineering Thermodynamics* — Ch. 3 (cubic EOS), Ch. 6 (residual properties), Ch. 10/13 (fugacity); App. B & C (property data).
