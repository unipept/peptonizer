import './style.css'
import javascriptLogo from './javascript.svg'
import peptonizerLogo from '/peptonizer.jpg'
import { peptonize } from './src/peptonizer.js'

document.querySelector('#app').innerHTML = `
  <div id="app">
    <div>
      <img src="${peptonizerLogo}" class="logo" alt="Peptonizer logo" />
      <a href="https://developer.mozilla.org/en-US/docs/Web/JavaScript" target="_blank">
        <img src="${javascriptLogo}" class="logo vanilla" alt="JavaScript logo" />
      </a>
      <h1>Peptonizer2000</h1>
      
      <div id="inputs">
        <div style="margin: 16px;">
            <span>1. Upload a PSM file</span>
            <label for="file-input" class="file-label">⇪ Choose PSM File</label>
            <input type="file" id="file-input" accept=".tsv" />
            <div id="file-input-label"></div>
        </div>
               
        <div style="margin: 8px;">
            <span>2. Start the Peptonizer!</span>    
            <button id="peptonize-button" disabled>↯ Start to Peptonize!</button>
        </div>
      </div>
      
      <div id="loading-spinner" hidden>
          <div class="lds-roller"><div></div><div></div><div></div><div></div><div></div><div></div><div></div><div></div></div>
          <div>Processing...</div>
          <div id="grid-progress-view"></div>
          <div id="graph-progress-view"></div>
          <div id="iterations-progress-view"></div>
          <div id="residual-progress-view"></div>
      </div>
      <div id="result-view" hidden>
          <h2>Final output</h2>
          <div id="peptonizer-chart" style="width:100%; min-width: 1000px; height:400px;"></div>
      </div>
    </div>
  </div>
`

let fileContents = "";

document.getElementById('file-input').addEventListener('change', (event) => {
    const file = event.target.files[0];
    if (file) {
        const reader = new FileReader();
        reader.onload = function(e) {
            fileContents = e.target.result;
            document.getElementById('peptonize-button').disabled = false; // Enable the button once the file is read
        }
        reader.readAsText(file);
        document.getElementById("file-input-label").innerHTML = "1 file selected"
    }
});


class ProgressListener {
    constructor(
        gridProgressEl,
        graphProgressEl,
        iterationsProgressEl,
        residualProgressEl
    ) {
        this.gridProgressEl = gridProgressEl;
        this.graphProgressEl = graphProgressEl;
        this.iterationsProgressEl = iterationsProgressEl;
        this.residualProgressEl = residualProgressEl;
    }

    gridUpdated(currentAlpha, currentBeta, currentPrior) {
        this.gridProgressEl.innerHTML = `Currently training graph with parameters α = ${currentAlpha}, β = ${currentBeta}, γ = ${currentPrior}.`
    }

    graphsUpdated(currentGraph, totalGraphs) {
        this.graphProgressEl.innerHTML = `Finished processing graph ${currentGraph} / ${totalGraphs}`;
    }

    maxResidualUpdated(maxResidual, tolerance) {
        this.residualProgressEl.innerHTML = `Improved maximum residual metric to ${maxResidual}. Tolerance is ${tolerance}`;
    }

    iterationsUpdated(currentIteration, totalIterations) {
        this.iterationsProgressEl.innerHTML = `Finished iteration ${currentIteration} / ${totalIterations}.`;
    }
}


const startToPeptonize = async function() {
    document.getElementById("result-view").hidden = true;
    document.getElementById("inputs").hidden = true;
    document.getElementById("loading-spinner").hidden = false;

    const alphas = [0.2, 0.5, 0.8];
    const betas = [0.2, 0.5, 0.8];
    const priors = [0.2, 0.5];

    const peptonizerResult = await peptonize(
        fileContents,
        alphas,
        betas,
        priors,
        new ProgressListener(
            document.getElementById("graph-progress-view"),
            document.getElementById("iterations-progress-view"),
            document.getElementById("residual-progress-view"),
        )
    );

    console.log(peptonizerResult);
    const entries = Object.entries(peptonizerResult[0]).map(([key, value]) => [key, parseFloat(value.toFixed(2))]);
    const sortedEntries = entries.sort((a, b) => b[1] - a[1]);


    // Extract keys and values from peptonizerResult
    const labels = sortedEntries.map(entry => entry[0]); // Sorted keys
    const values = sortedEntries.map(entry => entry[1]);

    // Render the chart with Highcharts
    Highcharts.chart('peptonizer-chart', {
        chart: {
            type: 'bar'
        },
        title: {
            text: 'Peptonizer Confidence Scores'
        },
        xAxis: {
            categories: labels.slice(0, 20),
            title: {
                text: 'Peptide IDs'
            }
        },
        yAxis: {
            min: 0,
            max: 1,
            title: {
                text: 'Confidence Score',
                align: 'high'
            },
            labels: {
                overflow: 'justify',
                format: '{value:.3f}'
            }
        },
        tooltip: {
            pointFormat: 'Confidence: <b>{point.y:.2f}</b>'
        },
        plotOptions: {
            bar: {
                dataLabels: {
                    enabled: true,
                    format: '{y:.3f}'
                }
            }
        },
        series: [{
            name: 'Confidence score',
            data: values.slice(0, 20)
        }]
    });

    document.getElementById("result-view").hidden = false;
    document.getElementById("loading-spinner").hidden = true;
}

document.getElementById("peptonize-button").addEventListener("click", startToPeptonize);
