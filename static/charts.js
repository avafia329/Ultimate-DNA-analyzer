function renderCountChart(data) {
    new Chart(document.getElementById('countChart'), {
        type: 'bar',
        data: {
            labels: Object.keys(data),
            datasets: [{
                label: 'Nucleotide Count',
                data: Object.values(data),
                backgroundColor: ['#ff9999', '#99ff99', '#9999ff', '#ffff99']
            }]
        }
    });
}

function renderGCChart(gcData) {
    new Chart(document.getElementById('gcChart'), {
        type: 'line',
        data: {
            labels: gcData.positions,
            datasets: [{
                label: 'GC Skew',
                data: gcData.values,
                borderColor: '#33cc33',
                fill: false
            }]
        }
    });
}
