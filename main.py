#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Flask Application for DUT Engine Visualization
"""

from flask import Flask, render_template_string, jsonify, request
import json
import sys
import os
from pathlib import Path
import threading
import time

# Importa o módulo DUT
import dut_cmb_scientific_engine as dut

app = Flask(__name__)

# Cache para resultados
results_cache = {
    'status': 'idle',
    'data': None,
    'error': None,
    'progress': 0
}

HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>DUT Scientific Engine 3.0 - NASA/ESA Grade</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }

        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #0f0c29, #302b63, #24243e);
            color: #ffffff;
            min-height: 100vh;
            padding: 20px;
        }

        .container {
            max-width: 1400px;
            margin: 0 auto;
        }

        header {
            text-align: center;
            padding: 40px 20px;
            background: rgba(255, 255, 255, 0.05);
            border-radius: 15px;
            margin-bottom: 30px;
            backdrop-filter: blur(10px);
        }

        h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
            background: linear-gradient(45deg, #00d4ff, #0099ff);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            background-clip: text;
        }

        .subtitle {
            color: #a0aec0;
            font-size: 1.1em;
        }

        .controls {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }

        .control-card {
            background: rgba(255, 255, 255, 0.08);
            padding: 20px;
            border-radius: 12px;
            backdrop-filter: blur(10px);
            border: 1px solid rgba(255, 255, 255, 0.1);
        }

        .control-card h3 {
            margin-bottom: 15px;
            color: #00d4ff;
            font-size: 1.1em;
        }

        button {
            width: 100%;
            padding: 12px 24px;
            background: linear-gradient(135deg, #00d4ff, #0099ff);
            color: white;
            border: none;
            border-radius: 8px;
            font-size: 1em;
            cursor: pointer;
            transition: all 0.3s ease;
            font-weight: 600;
        }

        button:hover {
            transform: translateY(-2px);
            box-shadow: 0 5px 20px rgba(0, 212, 255, 0.4);
        }

        button:disabled {
            background: #4a5568;
            cursor: not-allowed;
            transform: none;
        }

        .loading {
            text-align: center;
            padding: 40px;
            background: rgba(255, 255, 255, 0.05);
            border-radius: 12px;
            margin: 20px 0;
        }

        .spinner {
            border: 4px solid rgba(255, 255, 255, 0.1);
            border-radius: 50%;
            border-top: 4px solid #00d4ff;
            width: 50px;
            height: 50px;
            animation: spin 1s linear infinite;
            margin: 0 auto 20px;
        }

        @keyframes spin {
            0% { transform: rotate(0deg); }
            100% { transform: rotate(360deg); }
        }

        .results-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }

        .result-card {
            background: rgba(255, 255, 255, 0.08);
            padding: 25px;
            border-radius: 12px;
            backdrop-filter: blur(10px);
            border: 1px solid rgba(255, 255, 255, 0.1);
        }

        .result-card h3 {
            color: #00d4ff;
            margin-bottom: 15px;
            font-size: 1.3em;
        }

        .metric {
            display: flex;
            justify-content: space-between;
            padding: 10px 0;
            border-bottom: 1px solid rgba(255, 255, 255, 0.1);
        }

        .metric:last-child {
            border-bottom: none;
        }

        .metric-label {
            color: #a0aec0;
        }

        .metric-value {
            font-weight: 600;
            color: #ffffff;
        }

        .metric-value.positive {
            color: #48bb78;
        }

        .metric-value.negative {
            color: #f56565;
        }

        .chart-container {
            background: rgba(255, 255, 255, 0.08);
            padding: 25px;
            border-radius: 12px;
            backdrop-filter: blur(10px);
            border: 1px solid rgba(255, 255, 255, 0.1);
            margin-bottom: 20px;
        }

        .chart-container h3 {
            color: #00d4ff;
            margin-bottom: 20px;
            font-size: 1.3em;
        }

        canvas {
            max-height: 400px;
        }

        .error {
            background: rgba(245, 101, 101, 0.2);
            border: 1px solid #f56565;
            padding: 20px;
            border-radius: 12px;
            margin: 20px 0;
        }

        .success {
            background: rgba(72, 187, 120, 0.2);
            border: 1px solid #48bb78;
            padding: 20px;
            border-radius: 12px;
            margin: 20px 0;
        }

        pre {
            background: rgba(0, 0, 0, 0.3);
            padding: 15px;
            border-radius: 8px;
            overflow-x: auto;
            font-size: 0.9em;
            line-height: 1.5;
        }

        .badge {
            display: inline-block;
            padding: 4px 12px;
            border-radius: 12px;
            font-size: 0.85em;
            font-weight: 600;
            margin-left: 10px;
        }

        .badge.success {
            background: #48bb78;
            color: white;
        }

        .badge.warning {
            background: #ed8936;
            color: white;
        }

        .badge.info {
            background: #4299e1;
            color: white;
        }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>🌌 DUT Scientific Engine 3.0</h1>
            <p class="subtitle">Dead Universe Theory • NASA/ESA Grade Analysis</p>
        </header>

        <div class="controls">
            <div class="control-card">
                <h3>Validation Mode</h3>
                <button onclick="runAnalysis('validate')">Run Validation</button>
            </div>

            <div class="control-card">
                <h3>Comparison Mode</h3>
                <button onclick="runAnalysis('compare')">DUT vs ΛCDM</button>
            </div>

            <div class="control-card">
                <h3>DESI BAO Mode</h3>
                <button onclick="runAnalysis('desi')">Run with DESI Data</button>
            </div>

            <div class="control-card">
                <h3>Dimensional Test</h3>
                <button onclick="runAnalysis('test')">Run Tests</button>
            </div>
        </div>

        <div id="loading" class="loading" style="display: none;">
            <div class="spinner"></div>
            <p>Processing cosmological data... This may take a few moments.</p>
            <p id="progress-text">Progress: 0%</p>
        </div>

        <div id="results"></div>
    </div>

    <script>
        let pollInterval = null;

        async function runAnalysis(mode) {
            const resultsDiv = document.getElementById('results');
            const loadingDiv = document.getElementById('loading');

            resultsDiv.innerHTML = '';
            loadingDiv.style.display = 'block';

            // Desabilita todos os botões
            document.querySelectorAll('button').forEach(btn => btn.disabled = true);

            try {
                const response = await fetch('/run', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ mode: mode })
                });

                const data = await response.json();

                if (data.status === 'started') {
                    // Inicia polling para verificar progresso
                    startPolling();
                } else {
                    displayResults(data);
                    loadingDiv.style.display = 'none';
                }
            } catch (error) {
                resultsDiv.innerHTML = `<div class="error">Error: ${error.message}</div>`;
                loadingDiv.style.display = 'none';
            } finally {
                document.querySelectorAll('button').forEach(btn => btn.disabled = false);
            }
        }

        function startPolling() {
            if (pollInterval) clearInterval(pollInterval);

            pollInterval = setInterval(async () => {
                try {
                    const response = await fetch('/status');
                    const data = await response.json();

                    document.getElementById('progress-text').textContent =
                        `Progress: ${data.progress}%`;

                    if (data.status === 'completed') {
                        clearInterval(pollInterval);
                        displayResults(data.data);
                        document.getElementById('loading').style.display = 'none';
                        document.querySelectorAll('button').forEach(btn => btn.disabled = false);
                    } else if (data.status === 'error') {
                        clearInterval(pollInterval);
                        document.getElementById('results').innerHTML =
                            `<div class="error">Error: ${data.error}</div>`;
                        document.getElementById('loading').style.display = 'none';
                        document.querySelectorAll('button').forEach(btn => btn.disabled = false);
                    }
                } catch (error) {
                    console.error('Polling error:', error);
                }
            }, 1000);
        }

        function displayResults(data) {
            const resultsDiv = document.getElementById('results');

            if (!data) {
                resultsDiv.innerHTML = '<div class="error">No data received</div>';
                return;
            }

            let html = '';

            // Validation Results
            if (data.validation) {
                html += generateValidationCard(data.validation);
            }

            // Test Results
            if (data.test) {
                html += generateTestCard(data);
            }

            // Comparison Results
            if (data.baseline && data.dut) {
                html += generateComparisonCards(data);
                html += generateCharts(data);
            }

            // DUT Validation
            if (data.dut_validation) {
                html += generateDUTValidationCard(data.dut_validation);
            }

            // Metadata
            if (data.computation_metadata) {
                html += generateMetadataCard(data.computation_metadata);
            }

            resultsDiv.innerHTML = html;
        }

        function generateValidationCard(validation) {
            const isValid = validation.parameters_valid;
            const statusClass = isValid ? 'success' : 'error';
            const statusText = isValid ? 'VALID' : 'INVALID';

            return `
                <div class="result-card">
                    <h3>Physical Validation <span class="badge ${isValid ? 'success' : 'warning'}">${statusText}</span></h3>
                    <div class="metric">
                        <span class="metric-label">G_eff/G Today:</span>
                        <span class="metric-value">${validation.geff_today?.toFixed(6) || 'N/A'}</span>
                    </div>
                    ${validation.warnings && validation.warnings.length > 0 ?
                        `<div style="margin-top: 15px;">
                            <strong>Warnings:</strong>
                            <pre>${validation.warnings.join('\\n')}</pre>
                        </div>` : ''}
                </div>
            `;
        }

        function generateTestCard(data) {
            const isPassed = data.test === 'PASSED';
            return `
                <div class="result-card">
                    <h3>Dimensional Consistency Test
                        <span class="badge ${isPassed ? 'success' : 'warning'}">
                            ${data.test}
                        </span>
                    </h3>
                    ${data.error ? `<pre>${data.error}</pre>` : ''}
                </div>
            `;
        }

        function generateComparisonCards(data) {
            const deltaχ2 = data.delta_chi2;
            const improvement = deltaχ2 < 0;

            return `
                <div class="results-grid">
                    <div class="result-card">
                        <h3>ΛCDM Baseline</h3>
                        <div class="metric">
                            <span class="metric-label">χ² Total:</span>
                            <span class="metric-value">${data.baseline.chi2?.toFixed(2) || 'N/A'}</span>
                        </div>
                    </div>

                    <div class="result-card">
                        <h3>DUT Model</h3>
                        <div class="metric">
                            <span class="metric-label">χ² Total:</span>
                            <span class="metric-value">${data.dut.chi2?.toFixed(2) || 'N/A'}</span>
                        </div>
                    </div>

                    <div class="result-card">
                        <h3>Model Comparison</h3>
                        <div class="metric">
                            <span class="metric-label">Δχ²:</span>
                            <span class="metric-value ${improvement ? 'positive' : 'negative'}">
                                ${deltaχ2?.toFixed(2) || 'N/A'}
                            </span>
                        </div>
                        <div class="metric">
                            <span class="metric-label">Improvement:</span>
                            <span class="metric-value ${improvement ? 'positive' : 'negative'}">
                                ${improvement ? 'YES ✓' : 'NO ✗'}
                            </span>
                        </div>
                    </div>
                </div>
            `;
        }

        function generateDUTValidationCard(validation) {
            return `
                <div class="result-card">
                    <h3>DUT Physics Validation</h3>
                    <div class="metric">
                        <span class="metric-label">φ Standard Deviation:</span>
                        <span class="metric-value">${validation.phi_std?.toFixed(6) || 'N/A'}</span>
                    </div>
                    <div class="metric">
                        <span class="metric-label">φ Final:</span>
                        <span class="metric-value">${validation.phi_final?.toFixed(6) || 'N/A'}</span>
                    </div>
                    <div class="metric">
                        <span class="metric-label">G_eff Final:</span>
                        <span class="metric-value">${validation.G_eff_final?.toFixed(6) || 'N/A'}</span>
                    </div>
                    <div class="metric">
                        <span class="metric-label">G_eff at z=10:</span>
                        <span class="metric-value">${validation.G_eff_z10?.toFixed(6) || 'N/A'}</span>
                    </div>
                    <div class="metric">
                        <span class="metric-label">G_eff Change (z=10→0):</span>
                        <span class="metric-value">${(validation.G_eff_fractional_change_z10_to_0 * 100)?.toFixed(3) || 'N/A'}%</span>
                    </div>
                </div>
            `;
        }

        function generateMetadataCard(metadata) {
            const modules = metadata.modules_enabled || {};
            const modulesList = Object.entries(modules)
                .map(([key, val]) => `<span class="badge ${val ? 'success' : 'warning'}">${key}</span>`)
                .join(' ');

            return `
                <div class="result-card">
                    <h3>Computation Metadata</h3>
                    <div class="metric">
                        <span class="metric-label">Engine Version:</span>
                        <span class="metric-value">${metadata.engine_version || 'N/A'}</span>
                    </div>
                    <div class="metric">
                        <span class="metric-label">Timestamp:</span>
                        <span class="metric-value">${metadata.utc_time || 'N/A'}</span>
                    </div>
                    ${metadata.dimensional_test ? `
                    <div class="metric">
                        <span class="metric-label">Dimensional Test:</span>
                        <span class="metric-value">${metadata.dimensional_test}</span>
                    </div>` : ''}
                    <div style="margin-top: 15px;">
                        <strong>Enabled Modules:</strong><br><br>
                        ${modulesList}
                    </div>
                </div>
            `;
        }

        function generateCharts(data) {
            if (!data.baseline || !data.dut) return '';

            return `
                <div class="chart-container">
                    <h3>Chi-Squared Comparison</h3>
                    <canvas id="chi2Chart"></canvas>
                </div>
                <script>
                    const ctx = document.getElementById('chi2Chart');
                    new Chart(ctx, {
                        type: 'bar',
                        data: {
                            labels: ['ΛCDM', 'DUT', 'Improvement (Δχ²)'],
                            datasets: [{
                                label: 'χ² Value',
                                data: [
                                    ${data.baseline.chi2 || 0},
                                    ${data.dut.chi2 || 0},
                                    Math.abs(${data.delta_chi2 || 0})
                                ],
                                backgroundColor: [
                                    'rgba(255, 99, 132, 0.7)',
                                    'rgba(54, 162, 235, 0.7)',
                                    'rgba(75, 192, 192, 0.7)'
                                ],
                                borderColor: [
                                    'rgba(255, 99, 132, 1)',
                                    'rgba(54, 162, 235, 1)',
                                    'rgba(75, 192, 192, 1)'
                                ],
                                borderWidth: 2
                            }]
                        },
                        options: {
                            responsive: true,
                            maintainAspectRatio: true,
                            scales: {
                                y: {
                                    beginAtZero: true,
                                    grid: {
                                        color: 'rgba(255, 255, 255, 0.1)'
                                    },
                                    ticks: {
                                        color: '#ffffff'
                                    }
                                },
                                x: {
                                    grid: {
                                        color: 'rgba(255, 255, 255, 0.1)'
                                    },
                                    ticks: {
                                        color: '#ffffff'
                                    }
                                }
                            },
                            plugins: {
                                legend: {
                                    labels: {
                                        color: '#ffffff'
                                    }
                                }
                            }
                        }
                    });

            `;
        }
    </script>
</body>
</html>
"""

def run_dut_analysis(mode='compare', use_desi=False):
    """Executa a análise DUT em uma thread separada"""
    try:
        results_cache['status'] = 'running'
        results_cache['progress'] = 10

        analysis = dut.DUTCompleteAnalysis()
        results_cache['progress'] = 30

        if use_desi:
            analysis.bao_data = dut.load_desi_bao_placeholder()

        results_cache['progress'] = 50

        if mode == 'validate':
            params = analysis.create_dut_parameters()
            warnings = params.validate_physics()
            results = {
                "validation": {
                    "parameters_valid": len(warnings) == 0,
                    "warnings": warnings,
                    "geff_today": 1.0 / (1.0 + params.xi * params.phi_ini**2)
                },
                "engine_version": "3.0"
            }
        elif mode == 'test':
            dut.test_dimensional_consistency()
            results = {"test": "PASSED", "engine_version": "3.0"}
        else:
            results = analysis.run_comparison()

        results_cache['progress'] = 100
        results_cache['status'] = 'completed'
        results_cache['data'] = results
        results_cache['error'] = None

    except Exception as e:
        results_cache['status'] = 'error'
        results_cache['error'] = str(e)
        results_cache['data'] = None

@app.route('/')
def index():
    return render_template_string(HTML_TEMPLATE)

@app.route('/run', methods=['POST'])
def run_analysis():
    data = request.get_json()
    mode = data.get('mode', 'compare')
    use_desi = mode == 'desi'

    # Inicia a análise em uma thread separada
    thread = threading.Thread(target=run_dut_analysis, args=(mode, use_desi))
    thread.daemon = True
    thread.start()

    return jsonify({'status': 'started'})

@app.route('/status')
def get_status():
    return jsonify({
        'status': results_cache['status'],
        'progress': results_cache['progress'],
        'data': results_cache['data'],
        'error': results_cache['error']
    })

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 8000))
    app.run(host='0.0.0.0', port=port, debug=False)