document.getElementById('dnaForm').addEventListener('submit', function(e) {
    e.preventDefault();
    const seq = document.getElementById('dnaSequence').value;
    fetch('/analyze', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ sequence: seq })
    })
    .then(res => res.json())
    .then(data => {
        document.getElementById('analysisText').textContent = data.text;
        renderCountChart(data.counts);
        renderGCChart(data.gc_skew);
    });
});

function fetchFeature(feature) {
    const seq = document.getElementById('dnaSequence').value;
    fetch(`/feature/${feature}`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ sequence: seq })
    })
    .then(res => res.json())
    .then(data => {
        document.getElementById('featureOutput').textContent = data.result;
    });
}

function exportPDF() {
    fetch('/export', { method: 'GET' })
    .then(res => res.blob())
    .then(blob => {
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = "dna_report.pdf";
        a.click();
    });
}
