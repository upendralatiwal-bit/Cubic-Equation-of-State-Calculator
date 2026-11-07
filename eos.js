// Generic cubic EOS engine (Smith, Van Ness & Abbott form)
// P = RT/(V-b) - a(T)/[(V+eps*b)(V+sig*b)]
const R = 8.314; // J/mol/K

const COMPONENTS = {
  // Tc [K], Pc [bar], omega, Cp/R = A + B T + C T^2 + D T^-2  (SVN Table C.1)
  "Methane":   {Tc:190.6, Pc:45.99, w:0.012, A:1.702, B:9.081e-3, C:-2.164e-6, D:0},
  "Ethane":    {Tc:305.3, Pc:48.72, w:0.100, A:1.131, B:19.225e-3, C:-5.561e-6, D:0},
  "Propane":   {Tc:369.8, Pc:42.48, w:0.152, A:1.213, B:28.785e-3, C:-8.824e-6, D:0},
  "n-Butane":  {Tc:425.1, Pc:37.96, w:0.200, A:1.935, B:36.915e-3, C:-11.402e-6, D:0},
  "Ethylene":  {Tc:282.3, Pc:50.40, w:0.087, A:1.424, B:14.394e-3, C:-4.392e-6, D:0},
  "Carbon dioxide": {Tc:304.2, Pc:73.83, w:0.224, A:5.457, B:1.045e-3, C:0, D:-1.157e5},
  "Carbon monoxide":{Tc:132.9, Pc:34.99, w:0.048, A:3.376, B:0.557e-3, C:0, D:-0.031e5},
  "Nitrogen":  {Tc:126.2, Pc:34.00, w:0.038, A:3.280, B:0.593e-3, C:0, D:0.040e5},
  "Oxygen":    {Tc:154.6, Pc:50.43, w:0.022, A:3.639, B:0.506e-3, C:0, D:-0.227e5},
  "Hydrogen":  {Tc:33.19, Pc:13.13, w:-0.216, A:3.249, B:0.422e-3, C:0, D:0.083e5},
  "Ammonia":   {Tc:405.7, Pc:112.80, w:0.253, A:3.578, B:3.020e-3, C:0, D:-0.186e5},
  "Water":     {Tc:647.1, Pc:220.55, w:0.345, A:3.470, B:1.450e-3, C:0, D:0.121e5}
};

const EOS = {
  vdW: {name:"van der Waals",     sig:0, eps:0, Om:1/8,     Psi:27/64},
  RK:  {name:"Redlich–Kwong",     sig:1, eps:0, Om:0.08664, Psi:0.42748},
  SRK: {name:"Soave–Redlich–Kwong", sig:1, eps:0, Om:0.08664, Psi:0.42748},
  PR:  {name:"Peng–Robinson",     sig:1+Math.SQRT2, eps:1-Math.SQRT2, Om:0.07780, Psi:0.45724}
};

// alpha(Tr) and dln(alpha)/dln(T)
function alphaFn(model, Tr, w) {
  if (model === "vdW") return {a:1, d:0};
  if (model === "RK")  return {a:Math.pow(Tr,-0.5), d:-0.5};
  const m = model === "SRK" ? 0.480 + 1.574*w - 0.176*w*w
                            : 0.37464 + 1.54226*w - 0.26992*w*w;
  const s = 1 + m*(1 - Math.sqrt(Tr));
  return {a:s*s, d:-m*Math.sqrt(Tr)/s};
}

// Real roots of z^3 + c2 z^2 + c1 z + c0 = 0
function cubicRoots(c2, c1, c0) {
  const p = c1 - c2*c2/3, q = 2*c2*c2*c2/27 - c2*c1/3 + c0;
  const disc = q*q/4 + p*p*p/27, sh = -c2/3;
  let r;
  if (disc > 1e-14) {
    const u = Math.cbrt(-q/2 + Math.sqrt(disc)), v = Math.cbrt(-q/2 - Math.sqrt(disc));
    r = [u + v + sh];
  } else if (Math.abs(disc) <= 1e-14) {
    const u = Math.cbrt(-q/2); r = [2*u + sh, -u + sh];
  } else {
    const rho = Math.sqrt(-p*p*p/27), th = Math.acos(Math.max(-1, Math.min(1, -q/(2*rho))));
    const m = 2*Math.cbrt(rho);
    r = [0,1,2].map(k => m*Math.cos((th + 2*Math.PI*k)/3) + sh);
  }
  return r.sort((a,b)=>a-b);
}

// Ideal-gas Cp integrals from T0 to T (Cp/R form)
function igH(c, T0, T) { // returns ΔH/R  [K]
  return c.A*(T-T0) + c.B/2*(T*T-T0*T0) + c.C/3*(T**3-T0**3) - c.D*(1/T - 1/T0);
}
function igS(c, T0, T) { // returns ΔS/R from Cp only
  return c.A*Math.log(T/T0) + c.B*(T-T0) + c.C/2*(T*T-T0*T0) - c.D/2*(1/(T*T) - 1/(T0*T0));
}

/**
 * Solve the EOS at T [K], P [bar] for component props c.
 * Returns every physical root with Z, V, HR, SR, phi, f, plus the stable root.
 */
function solve(model, c, T, P, ref = {T0:298.15, P0:1}) {
  const e = EOS[model];
  const Tr = T/c.Tc, Pr = P/c.Pc;
  const al = alphaFn(model, Tr, c.w);
  const beta = e.Om*Pr/Tr;                   // = bP/RT
  const q = e.Psi*al.a/(e.Om*Tr);            // = a/(bRT)
  const Aa = q*beta;                         // = aP/(RT)^2
  const B = beta, s = e.sig, ep = e.eps;
  const c2 = (s+ep)*B - 1 - B;
  const c1 = Aa + ep*s*B*B - (s+ep)*B*(1+B);
  const c0 = -(Aa*B + ep*s*B*B*(1+B));
  const Zs = cubicRoots(c2, c1, c0).filter(z => z > B + 1e-12);

  const a_SI = e.Psi*al.a*R*R*c.Tc*c.Tc/(c.Pc*1e5);  // Pa m6/mol2
  const b_SI = e.Om*R*c.Tc/(c.Pc*1e5);               // m3/mol

  const dHig = R*igH(c, ref.T0, T);
  const dSig = R*(igS(c, ref.T0, T) - Math.log(P/ref.P0));

  const roots = Zs.map(Z => {
    const I = (s === ep) ? B/(Z + ep*B) : Math.log((Z + s*B)/(Z + ep*B))/(s - ep);
    const HR = R*T*(Z - 1 + (al.d - 1)*q*I);      // J/mol
    const SR = R*(Math.log(Z - B) + al.d*q*I);    // J/mol/K
    const lnphi = Z - 1 - Math.log(Z - B) - q*I;
    const phi = Math.exp(lnphi);
    return {Z, V: Z*R*T/(P*1e5)*1e6 /* cm3/mol */, HR, SR, phi, f: phi*P,
            H: dHig + HR, S: dSig + SR, G_R: lnphi};
  });
  // Stable root = lowest Gibbs energy = lowest ln(phi)
  let stable = roots.reduce((m, r) => (r.G_R < m.G_R ? r : m), roots[0]);
  return {Tr, Pr, alpha:al.a, a:a_SI, b:b_SI, A:Aa, B, roots, stable, dHig, dSig,
          phase: roots.length > 1 ? "two-phase region (3 roots)" :
                 (T > c.Tc ? "supercritical / gas" : (stable.Z > 0.3 ? "vapour" : "liquid"))};
}

// Pressure from EOS at T [K], V [cm3/mol] -> bar
function pressure(model, c, T, Vcm3) {
  const e = EOS[model], al = alphaFn(model, T/c.Tc, c.w);
  const a = e.Psi*al.a*R*R*c.Tc*c.Tc/(c.Pc*1e5), b = e.Om*R*c.Tc/(c.Pc*1e5);
  const V = Vcm3*1e-6;
  return (R*T/(V - b) - a/((V + e.eps*b)*(V + e.sig*b)))/1e5;
}

if (typeof module !== "undefined") module.exports = {solve, pressure, COMPONENTS, EOS};
