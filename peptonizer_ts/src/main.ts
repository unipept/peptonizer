import './style.css'
import typescriptLogo from "./typescript.svg"
import peptonizerLogo from "./peptonizer.jpg"
import { GridSearchProgressListener } from "./GridSearchProgressListener.ts";
import { Peptonizer } from "./Peptonizer.ts";

document.querySelector<HTMLDivElement>('#app')!.innerHTML= `
  <div id="app">
    <div>
      <img src="${peptonizerLogo}" class="logo" alt="Peptonizer logo" />
      <a href="https://developer.mozilla.org/en-US/docs/Web/JavaScript" target="_blank">
        <img src="${typescriptLogo}" class="logo vanilla" alt="JavaScript logo" />
      </a>
      <h1>Peptonizer2000</h1>
      
      <div id="inputs">
        <div style="margin: 16px;">
            <span>1. Upload a PSM file</span>
            <label for="file-input" class="file-label">⇪ Choose PSM File</label>
            <input type="file" id="file-input" accept=".tsv,.txt" />
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
          <div id="progress-view" style="display: flex;"></div>
      </div>
      <div id="result-view" hidden>
          <h2>Final output</h2>
          <div id="peptonizer-chart" style="width:100%; min-width: 1000px; height:400px;"></div>
      </div>
    </div>
  </div>
`

let fileContents = "";

document.querySelector<HTMLDivElement>('#file-input')!.addEventListener('change', (event: Event) => {
    const input = event.target as HTMLInputElement;
    const file = input.files ? input.files[0] : null;

    if (file) {
        const reader = new FileReader();
        reader.onload = function(e: ProgressEvent<FileReader>) {
            fileContents = e.target?.result as string;
            document.querySelector<HTMLButtonElement>('#peptonize-button')!.disabled = false; // Enable the button once the file is read
        }
        reader.readAsText(file);
        document.getElementById("file-input-label")!.innerHTML = "1 file selected"
    }
});


type ProgressViewContainer = {
    gridProgressView: HTMLDivElement,
    graphProgressView: HTMLDivElement,
    residualProgressView: HTMLDivElement,
    iterationsProgressView: HTMLDivElement
}

class ProgressListener implements GridSearchProgressListener {
    private progressViews: ProgressViewContainer[] = [];

    constructor(
        private progressView: HTMLElement,
        private workerCount: number
    ) {
        this.progressViews = [];
        progressView.innerHTML = '';
        for (let worker = 0; worker < this.workerCount; worker++) {
            const div = document.createElement('div');
            div.className += 'worker-status';
            const workerView = document.createElement('div');
            workerView.innerHTML = `Status for worker ${worker}`;
            const gridProgressView = document.createElement('div');
            const graphProgressView = document.createElement('div');
            const residualProgressView = document.createElement('div');
            const iterationsProgressView = document.createElement('div');

            div.appendChild(workerView);
            div.appendChild(gridProgressView);
            div.appendChild(graphProgressView);
            div.appendChild(residualProgressView);
            div.appendChild(iterationsProgressView);

            this.progressView.appendChild(div);

            this.progressViews.push({
                gridProgressView,
                graphProgressView,
                residualProgressView,
                iterationsProgressView
            });
        }

    }

    gridUpdated(
        currentAlpha: number,
        currentBeta: number,
        currentPrior: number,
        workerId: number
    ) {
        this.progressViews[workerId].gridProgressView.innerHTML = `Currently training model with parameters α = ${currentAlpha}, β = ${currentBeta}, γ = ${currentPrior}.`
    }

    graphsUpdated(
        currentGraph: number,
        totalGraphs: number,
        workerId: number
    ) {
        this.progressViews[workerId].graphProgressView.innerHTML = `Finished processing graph ${currentGraph} / ${totalGraphs}`;
    }

    maxResidualUpdated(
        maxResidual: number,
        tolerance: number,
        workerId: number
    ) {
        this.progressViews[workerId].residualProgressView.innerHTML = `Improved maximum residual metric to ${maxResidual}. Tolerance is ${tolerance}`;
    }

    iterationsUpdated(
        currentIteration: number,
        totalIterations: number,
        workerId: number
    ) {
        this.progressViews[workerId].iterationsProgressView.innerHTML = `Finished iteration ${currentIteration} / ${totalIterations}.`;
    }
}


const startToPeptonize = async function() {
    const resultView: HTMLElement = document.getElementById("result-view")!;
    const inputElement: HTMLElement = document.getElementById("inputs")!;
    const loadingSpinner: HTMLElement = document.getElementById("loading-spinner")!;

    resultView.hidden = true;
    inputElement.hidden = true;
    loadingSpinner.hidden = false;

    const alphas = [0.2, 0.5, 0.8];
    const betas = [0.2, 0.5, 0.8];
    const priors = [0.2, 0.5];

    const peptonizer = new Peptonizer();

    const peptonizerResult = await peptonizer.peptonize(
        fileContents,
        alphas,
        betas,
        priors,
        new ProgressListener(document.getElementById("progress-view")!, 1)
    );

    console.log(peptonizerResult);
    const entries = Object.entries(peptonizerResult[0]).map(([key, value]) => [key, parseFloat(value.toFixed(2))]);
    // @ts-ignore
    const sortedEntries = entries.sort((a, b) => b[1] - a[1]);


    // Extract keys and values from peptonizerResult
    const labels = sortedEntries.map(entry => entry[0]); // Sorted keys
    const values = sortedEntries.map(entry => entry[1]);

    // Render the chart with Highcharts
    // @ts-ignore
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

    resultView.hidden = false;
    loadingSpinner.hidden = true;
    inputElement.hidden = false;
}

document.getElementById("peptonize-button")!.addEventListener("click", startToPeptonize);
