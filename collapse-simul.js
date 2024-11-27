// Define the parameters and initial conditions
const parameterDefinitions = [
    {name: 'alpha_C', label: 'Alpha_C', default: 1.0, tooltip: 'Weight for complexity in utility'},
    {name: 'alpha_S', label: 'Alpha_S', default: 1.0, tooltip: 'Weight for specialization in complexity'},
    {name: 'beta_G', label: 'Beta_G', default: 1.2, tooltip: 'Weight for governance in complexity'},
    {name: 'gamma_comp', label: 'Gamma_comp', default: 0.5, tooltip: 'Weight for compensatory entropy in complexity'},
    {name: 'gamma_comms', label: 'Gamma_comms', default: 0.4, tooltip: 'Weight for communication entropy in utility'},
    {name: 'gamma_H', label: 'Gamma_H', default: 1.0, tooltip: 'Weight for entropy in utility'},
    {name: 'gamma_compute', label: 'Gamma_compute', default: 0.3, tooltip: 'Weight for computational entropy in utility'},
    {name: 'lambda_', label: 'Lambda', default: 0.1, tooltip: 'Decay rate for specialization overlaps'},
    {name: 'mu', label: 'Mu', default: 0.01, tooltip: 'Growth rate for specialization overlaps'},
    {name: 'gamma', label: 'Gamma', default: 0.2, tooltip: 'Growth rate for governance weights'},
    {name: 'delta', label: 'Delta', default: 0.01, tooltip: 'Decay rate for governance weights'},
    {name: 'beta_J', label: 'Beta_J', default: 1.0, tooltip: 'Governance influence on entropy'},
    {name: 'eta_comms', label: 'Eta_comms', default: 0.05, tooltip: 'Entropy generation rate from communications'},
    {name: 'theta_comms', label: 'Theta_comms', default: 0.02, tooltip: 'Entropy mitigation rate for communications'},
    {name: 'eta_compute', label: 'Eta_compute', default: 0.03, tooltip: 'Entropy generation rate from computations'},
    {name: 'theta_compute', label: 'Theta_compute', default: 0.01, tooltip: 'Entropy mitigation rate for computations'},
    {name: 'eta_units', label: 'Eta_units', default: 0.04, tooltip: 'Entropy generation rate from units'},
    {name: 'theta_units', label: 'Theta_units', default: 0.015, tooltip: 'Entropy mitigation rate for units'},
    {name: 'sigma_S_0', label: 'Sigma_S_0', default: 0.01, tooltip: 'Base noise level for specialization'},
    {name: 'sigma_G_0', label: 'Sigma_G_0', default: 0.01, tooltip: 'Base noise level for governance'},
    {name: 'sigma_H_comms_0', label: 'Sigma_H_comms_0', default: 0.005, tooltip: 'Base noise level for communication entropy'},
    {name: 'sigma_H_0', label: 'Sigma_H_0', default: 0.005, tooltip: 'Base noise level for total entropy'},
    {name: 'sigma_H_compute_0', label: 'Sigma_H_compute_0', default: 0.005, tooltip: 'Base noise level for computational entropy'},
    {name: 'sigma_w_jk_0', label: 'Sigma_w_jk_0', default: 0.005, tooltip: 'Base noise level for overlaps'},
    {name: 'sigma_w_j0_0', label: 'Sigma_w_j0_0', default: 0.005, tooltip: 'Base noise level for governance weights'}
];

const initialConditionDefinitions = [
    {name: 'S0', label: 'Initial Specialization (S0)', default: 10.0},
    {name: 'G0', label: 'Initial Governance (G0)', default: 5.0},
    {name: 'H_comms0', label: 'Initial Communication Entropy (H_comms0)', default: 2.0},
    {name: 'H0', label: 'Initial Total Entropy (H0)', default: 5.0},
    {name: 'H_compute0', label: 'Initial Computational Entropy (H_compute0)', default: 3.0},
    {name: 'w_jk0', label: 'Initial Overlap (w_jk0)', default: 0.3},
    {name: 'w_j0_0', label: 'Initial Governance Weight (w_j0_0)', default: 0.5}
];

// Time span for simulation
let t_start = 0;
let t_end = 100;
let dt = 0.01;

// Generate parameter input elements
function generateParameterInputs() {
    const parametersDiv = document.getElementById('parameters-inputs');
    const initialConditionsDiv = document.getElementById('initial-conditions-inputs');

    parameterDefinitions.forEach(param => {
        const group = document.createElement('div');
        group.className = 'parameter-group';

        const label = document.createElement('label');
        label.setAttribute('for', param.name);
        label.className = 'form-label';
        label.innerHTML = `${param.label} <i class="fas fa-info-circle text-secondary" data-bs-toggle="tooltip" data-bs-placement="right" title="${param.tooltip}"></i>`;

        const input = document.createElement('input');
        input.type = 'number';
        input.step = 'any';
        input.className = 'form-control';
        input.id = param.name;
        input.value = param.default;
        input.required = true;

        group.appendChild(label);
        group.appendChild(input);
        parametersDiv.appendChild(group);
    });

    initialConditionDefinitions.forEach(cond => {
        const group = document.createElement('div');
        group.className = 'parameter-group';

        const label = document.createElement('label');
        label.setAttribute('for', cond.name);
        label.className = 'form-label';
        label.textContent = cond.label;

        const input = document.createElement('input');
        input.type = 'number';
        input.step = 'any';
        input.className = 'form-control';
        input.id = cond.name;
        input.value = cond.default;
        input.required = true;

        group.appendChild(label);
        group.appendChild(input);
        initialConditionsDiv.appendChild(group);
    });

    // Initialize Bootstrap tooltips
    const tooltipTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle="tooltip"]'))
    tooltipTriggerList.map(function (tooltipTriggerEl) {
        return new bootstrap.Tooltip(tooltipTriggerEl)
    });
}

// Generate standard normal random number using Box-Muller transform
function randn_bm() {
    let u = 0, v = 0;
    while(u === 0) u = Math.random(); // Convert [0,1) to (0,1)
    while(v === 0) v = Math.random();
    return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
}

// Main simulation function using Euler-Maruyama method
function eulerMaruyamaIntegration(params, initial_conditions, R, t) {
    const num_w_jk = (R * (R - 1)) / 2;
    const num_w_j0 = R;
    const num_states = 5 + num_w_jk + num_w_j0;
    let y = [];

    // Initialize state variables
    let y0 = new Array(num_states).fill(0);
    y0[0] = initial_conditions.S0;
    y0[1] = initial_conditions.G0;
    y0[2] = initial_conditions.H_comms0;
    y0[3] = initial_conditions.H0;
    y0[4] = initial_conditions.H_compute0;
    for(let i = 0; i < num_w_jk; i++) {
        y0[5 + i] = initial_conditions.w_jk0;
    }
    for(let i = 0; i < num_w_j0; i++) {
        y0[5 + num_w_jk + i] = initial_conditions.w_j0_0;
    }
    y.push(y0);

    const k_sigma = 0.001; // Rate of volatility increase
    let collapseTime = null;

    for(let i = 1; i < t.length; i++) {
        const current_y = y[i - 1].slice();

        const S = current_y[0];
        const G = current_y[1];
        const H_comms = current_y[2];
        const H = current_y[3];
        const H_compute = current_y[4];
        const w_jk = current_y.slice(5, 5 + num_w_jk);
        const w_j0 = current_y.slice(5 + num_w_jk);

        // Governance Dynamics
        let sum_gw0 = 0;
        for(let j = 0; j < num_w_j0; j++) {
            sum_gw0 += params.gamma * w_j0[j] - params.delta * Math.pow(w_j0[j], 2);
        }

        // Entropy Contributions
        const dH_comms_dt = params.eta_comms * G - params.theta_comms * H_comms;
        const average_overlap = w_jk.reduce((a, b) => a + b, 0) / w_jk.length;
        const dH_units_dt = params.eta_units * S - params.theta_units * average_overlap;
        const dH_compute_dt = params.eta_compute * G - params.theta_compute * H_compute;
        const dH_dt = -params.beta_J * sum_gw0 + dH_comms_dt + dH_units_dt + dH_compute_dt;

        // Governance Growth
        const dG_dt = sum_gw0;

        // Specialization Dynamics
        const dS_dt = params.mu * (1 - average_overlap) - params.lambda_ * average_overlap;

        // Overlaps Between Jurisdictions Dynamics
        const dw_jk_dt = w_jk.map(w => -params.lambda_ * w + params.mu * (1 - w));

        // Governance Weights to Central Planner Dynamics
        const dw_j0_dt = w_j0.map(w => params.gamma * w - params.delta * Math.pow(w, 2));

        // Update state variables using Euler method
        let S_new = S + dS_dt * dt;
        let G_new = G + dG_dt * dt;
        let H_comms_new = H_comms + dH_comms_dt * dt;
        let H_new = H + dH_dt * dt;
        let H_compute_new = H_compute + dH_compute_dt * dt;
        let w_jk_new = w_jk.map((w, idx) => w + dw_jk_dt[idx] * dt);
        let w_j0_new = w_j0.map((w, idx) => w + dw_j0_dt[idx] * dt);

        // Update sigma values to be non-stationary (increasing)
        const sigma_S_t = params.sigma_S_0 * (1 + k_sigma * t[i]);
        const sigma_G_t = params.sigma_G_0 * (1 + k_sigma * t[i]);
        const sigma_H_comms_t = params.sigma_H_comms_0 * (1 + k_sigma * t[i]);
        const sigma_H_t = params.sigma_H_0 * (1 + k_sigma * t[i]);
        const sigma_H_compute_t = params.sigma_H_compute_0 * (1 + k_sigma * t[i]);
        const sigma_w_jk_t = params.sigma_w_jk_0 * (1 + k_sigma * t[i]);
        const sigma_w_j0_t = params.sigma_w_j0_0 * (1 + k_sigma * t[i]);

        // Sample noise with time-dependent (increasing) sigma
        const noise_S = randn_bm() * sigma_S_t * Math.sqrt(dt);
        const noise_G = randn_bm() * sigma_G_t * Math.sqrt(dt);
        const noise_H_comms = randn_bm() * sigma_H_comms_t * Math.sqrt(dt);
        const noise_H = randn_bm() * sigma_H_t * Math.sqrt(dt);
        const noise_H_compute = randn_bm() * sigma_H_compute_t * Math.sqrt(dt);
        const noise_w_jk = w_jk_new.map(() => randn_bm() * sigma_w_jk_t * Math.sqrt(dt));
        const noise_w_j0 = w_j0_new.map(() => randn_bm() * sigma_w_j0_t * Math.sqrt(dt));

        // Apply noise
        S_new += noise_S;
        G_new += noise_G;
        H_comms_new += noise_H_comms;
        H_new += noise_H;
        H_compute_new += noise_H_compute;
        w_jk_new = w_jk_new.map((w, idx) => w + noise_w_jk[idx]);
        w_j0_new = w_j0_new.map((w, idx) => w + noise_w_j0[idx]);

        // Prevent negative values
        S_new = Math.max(S_new, 0);
        G_new = Math.max(G_new, 0);
        H_comms_new = Math.max(H_comms_new, 0);
        H_new = Math.max(H_new, 0);
        H_compute_new = Math.max(H_compute_new, 0);
        w_jk_new = w_jk_new.map(w => Math.max(w, 0));
        w_j0_new = w_j0_new.map(w => Math.max(w, 0));

        // Calculate Coordination Benefits
        const coordination_benefits = params.lambda_ * w_jk_new.reduce((a, b) => a + b, 0) / 2;

        // Calculate Complexity
        const C_new = params.alpha_S * S_new + params.beta_G * G_new + params.gamma_comp * H_comms_new;

        // Calculate Utility
        const U_new = (params.alpha_C * C_new) +
                     (params.beta_G * G_new) +
                     coordination_benefits -
                     (params.gamma_comms * H_comms_new) -
                     (params.gamma_H * H_new) -
                     (params.gamma_compute * H_compute_new);

        // Update the state vector
        const y_new = [S_new, G_new, H_comms_new, H_new, H_compute_new, ...w_jk_new, ...w_j0_new];
        y.push(y_new);

        // Check for collapse condition
        if (U_new <= 0 && !collapseTime) {
            console.log(`Utility reached zero at time t = ${t[i].toFixed(2)}. Simulation halted.`);
            collapseTime = i;
            break;
        }
    }

    // Trim t if simulation halted early
    if (collapseTime !== null) {
        t = t.slice(0, collapseTime + 1);
        y = y.slice(0, collapseTime + 1);
    }

    return {y: y, t: t};
}

// Calculate Utility and Complexity
function calculateUtility(y, params, R) {
    const num_w_jk = (R * (R - 1)) / 2;
    const S = y.map(state => state[0]);
    const G = y.map(state => state[1]);
    const H_comms = y.map(state => state[2]);
    const H = y.map(state => state[3]);
    const H_compute = y.map(state => state[4]);
    const w_jk = y.map(state => state.slice(5, 5 + num_w_jk));

    const coordination_benefits = w_jk.map(wjks => params.lambda_ * wjks.reduce((a, b) => a + b, 0) / 2);
    const C = S.map((S_val, idx) => params.alpha_S * S_val + params.beta_G * G[idx] + params.gamma_comp * H_comms[idx]);
    const U = C.map((C_val, idx) => 
        (params.alpha_C * C_val) +
        (params.beta_G * G[idx]) +
        coordination_benefits[idx] -
        (params.gamma_comms * H_comms[idx]) -
        (params.gamma_H * H[idx]) -
        (params.gamma_compute * H_compute[idx])
    );

    return {U: U, C: C};
}

// Find Saddle Point
function findSaddlePoint(U, H) {
    const max_U_index = U.reduce((maxIdx, currentValue, currentIndex, array) => {
        return currentValue > array[maxIdx] ? currentIndex : maxIdx;
    }, 0);

    for(let i = max_U_index + 1; i < U.length; i++) {
        if(U[i] < U[i - 1]) {
            return {H_crit: H[i - 1], index: i - 1};
        }
    }
    return {H_crit: H[H.length - 1], index: H.length - 1};
}

// Plot Results using Plotly.js
function plotResults(t, y, U, C, R, H_crit, H_crit_index, params) {
    const num_w_jk = (R * (R - 1)) / 2;

    // Extract variables for plotting
    const S = y.map(state => state[0]);
    const G = y.map(state => state[1]);
    const H_comms = y.map(state => state[2]);
    const H = y.map(state => state[3]);
    const H_compute = y.map(state => state[4]);
    const w_jk = y.map(state => state.slice(5, 5 + num_w_jk));
    const w_j0 = y.map(state => state.slice(5 + num_w_jk));

    const coordination_benefits = w_jk.map(wjks => params.lambda_ * wjks.reduce((a, b) => a + b, 0) / 2);

    // Calculate Unit Entropy
    const H_units = H.map((H_val, idx) => H_val - H_comms[idx] - H_compute[idx]);

    // Calculate Governance Weights (Average) across Jurisdictions
    const governance_weights_avg = w_j0.map(weights => {
        const sum = weights.reduce((a, b) => a + b, 0);
        return sum / weights.length;
    });

    // Rename coordination_benefits to "Decentralized Coordination"
    const decentralized_coordination = coordination_benefits;

    // Clear previous plots
    const plotsDiv = document.getElementById('plots');
    plotsDiv.innerHTML = ''; // Clear previous plots

    // Create separate containers for each plot
    const plotIds = [
        'plot1', // Specialization, Governance, Complexity, Governance Weights (Average)
        'plot2', // Entropy Components, Decentralized Coordination
        'plot3', // Utility Over Time
        'plot4'  // Utility vs. Entropy (Phase Diagram)
    ];

    plotIds.forEach(id => {
        const plotContainer = document.createElement('div');
        plotContainer.id = id;
        plotContainer.style.width = '100%';
        plotContainer.style.height = '600px';
        plotContainer.className = 'mb-5';
        plotsDiv.appendChild(plotContainer);
    });

    // 1. Specialization, Governance, Complexity, and Governance Weights (Average) Over Time
    const traceS = { x: t, y: S, mode: 'lines', name: 'Specialization (S(t))', line: { color: '#1f77b4' } };
    const traceG = { x: t, y: G, mode: 'lines', name: 'Governance (G(t))', line: { color: '#ff7f0e' } };
    const traceC = { x: t, y: C, mode: 'lines', name: 'Complexity (C(t))', line: { color: '#2ca02c' } };
    const traceGW_avg = { x: t, y: governance_weights_avg, mode: 'lines', name: 'Governance Weights (Average)', line: { color: '#8c564b', dash: 'dash' } };
    const data1 = [traceS, traceG, traceC, traceGW_avg];
    const layout1 = {
        title: 'Specialization, Governance, Complexity, and Governance Weights Over Time',
        xaxis: { title: 'Time' },
        yaxis: { title: 'Value' },
        legend: { orientation: 'h', x: 0.5, xanchor: 'center' },
        margin: { t: 50 }
    };
    Plotly.newPlot('plot1', data1, layout1, {responsive: true});

    // 2. Entropy Components and Decentralized Coordination Over Time
    const traceH = { x: t, y: H, mode: 'lines', name: 'Total Entropy (H(t))', line: { color: '#d62728' } };
    const traceH_comms = { x: t, y: H_comms, mode: 'lines', name: 'Communication Entropy (H_comms(t))', line: { color: '#9467bd' } };
    const traceH_compute = { x: t, y: H_compute, mode: 'lines', name: 'Computational Entropy (H_compute(t))', line: { color: '#8c564b' } };
    const traceH_units = { x: t, y: H_units, mode: 'lines', name: 'Unit Entropy (H_units(t))', line: { color: '#e377c2', dash: 'dash' } };
    const traceDecentralizedCoordination = { x: t, y: decentralized_coordination, mode: 'lines', name: 'Decentralized Coordination', line: { color: '#7f7f7f', dash: 'dot' } };
    const traceH_crit = { 
        x: [t[0], t[t.length - 1]], 
        y: [H_crit, H_crit], 
        mode: 'lines', 
        name: 'Critical Entropy (H_crit)', 
        line: { dash: 'dashdot', color: '#000000' }
    };
    const data2 = [traceH, traceH_comms, traceH_compute, traceH_units, traceDecentralizedCoordination, traceH_crit];
    const layout2 = {
        title: 'Entropy Components and Decentralized Coordination Over Time',
        xaxis: { title: 'Time' },
        yaxis: { title: 'Entropy / Coordination' },
        legend: { orientation: 'h', x: 0.5, xanchor: 'center' },
        margin: { t: 50 }
    };
    Plotly.newPlot('plot2', data2, layout2, {responsive: true});

    // 3. Utility Over Time
    const traceU = { x: t, y: U, mode: 'lines', name: 'Utility (U(t))', line: { color: '#17becf' } };
    const traceZero = { 
        x: [t[0], t[t.length - 1]], 
        y: [0, 0], 
        mode: 'lines', 
        name: 'Zero Utility', 
        line: { dash: 'dash', color: '#000000' } 
    };
    const data3 = [traceU, traceZero];
    const layout3 = {
        title: 'Utility Over Time',
        xaxis: { title: 'Time' },
        yaxis: { title: 'Utility' },
        legend: { orientation: 'h', x: 0.5, xanchor: 'center' },
        margin: { t: 50 }
    };
    Plotly.newPlot('plot3', data3, layout3, {responsive: true});

    // 4. Utility vs. Entropy (Phase Diagram)
    const tracePhase = { x: H, y: U, mode: 'lines', name: 'Utility vs. Entropy', line: { color: '#bcbd22' } };
    const traceHcritPhase = { 
        x: [H_crit, H_crit], 
        y: [Math.min(...U), Math.max(...U)], 
        mode: 'lines', 
        name: 'Critical Entropy (H_crit)', 
        line: { dash: 'dashdot', color: '#7f7f7f' } 
    };
    const data4 = [tracePhase, traceHcritPhase];
    const layout4 = {
        title: 'Utility vs. Entropy (Phase Diagram)',
        xaxis: { title: 'Total Entropy (H(t))' },
        yaxis: { title: 'Utility (U(t))' },
        legend: { orientation: 'h', x: 0.5, xanchor: 'center' },
        margin: { t: 50 }
    };
    Plotly.newPlot('plot4', data4, layout4, {responsive: true});
}

// Display Alert Messages
function showAlert(message, type='info') {
    const alertContainer = document.getElementById('alert-container');
    const alert = document.createElement('div');
    alert.className = `alert alert-${type} alert-dismissible fade show`;
    alert.role = 'alert';
    alert.innerHTML = `
        ${message}
        <button type="button" class="btn-close" data-bs-dismiss="alert" aria-label="Close"></button>
    `;
    alertContainer.appendChild(alert);
    // Automatically remove the alert after 5 seconds
    setTimeout(() => {
        alert.classList.remove('show');
        alert.classList.add('hide');
    }, 5000);
}

// Main function to run the simulation
function runSimulation(event) {
    event.preventDefault();

    // Show loading overlay
    document.getElementById('loading-overlay').style.display = 'flex';

    // Read parameters from input elements
    let params = {};
    try {
        parameterDefinitions.forEach(param => {
            const value = parseFloat(document.getElementById(param.name).value);
            if (isNaN(value)) {
                throw new Error(`Invalid input for ${param.label}. Please enter a valid number.`);
            }
            params[param.name] = value;
        });

        let initial_conditions = {};
        initialConditionDefinitions.forEach(cond => {
            const value = parseFloat(document.getElementById(cond.name).value);
            if (isNaN(value)) {
                throw new Error(`Invalid input for ${cond.label}. Please enter a valid number.`);
            }
            initial_conditions[cond.name] = value;
        });

        // Read number of jurisdictions R
        const R_input = parseInt(document.getElementById('R').value);
        if (isNaN(R_input) || R_input < 2) {
            throw new Error('Number of Jurisdictions (R) must be an integer ≥ 2.');
        }
        const R = R_input;

        // Time array
        let t = [];
        for(let tt = t_start; tt <= t_end; tt += dt) {
            t.push(tt);
        }

        // Run the simulation
        const result = eulerMaruyamaIntegration(params, initial_conditions, R, t);
        let y = result.y;
        let t_sim = result.t;

        // Calculate Utility and Complexity
        const utilityResult = calculateUtility(y, params, R);
        const U = utilityResult.U;
        const C = utilityResult.C;

        // Determine H_crit using the saddle point after simulation
        const H = y.map(state => state[3]);
        const saddlePointResult = findSaddlePoint(U, H);
        const H_crit = saddlePointResult.H_crit;
        const H_crit_index = saddlePointResult.index;

        // Display critical entropy information
        if(saddlePointResult.index < U.length - 1) {
            showAlert(`Critical Entropy (H_crit) reached at time t = ${t_sim[H_crit_index].toFixed(2)} with H_crit = ${H_crit.toFixed(2)}.`, 'warning');
        } else {
            showAlert(`No organizational collapse detected within the simulation time span.`, 'success');
        }

        // Plot the results
        plotResults(t_sim, y, U, C, R, H_crit, H_crit_index, params);
    } catch (error) {
        showAlert(error.message, 'danger');
        console.error(error);
    } finally {
        // Hide loading overlay
        document.getElementById('loading-overlay').style.display = 'none';
    }
}

// Initialize the application
document.addEventListener('DOMContentLoaded', () => {
    generateParameterInputs();
    document.getElementById('parameters-form').addEventListener('submit', runSimulation);
});
