import './style.css'
import javascriptLogo from './javascript.svg'
import peptonizerLogo from '/peptonizer.jpg'
import { peptonize } from './src/peptonizer.js'

document.querySelector('#app').innerHTML = `
  <div>
    <img src="${peptonizerLogo}" class="logo" alt="Peptonizer logo" />
    <a href="https://developer.mozilla.org/en-US/docs/Web/JavaScript" target="_blank">
      <img src="${javascriptLogo}" class="logo vanilla" alt="JavaScript logo" />
    </a>
    <h1>Peptonizer2000</h1>
    <button id="peptonize-button">↯ Start to Peptonize!</button>
    <div id="loading-spinner" hidden>
        <div class="lds-roller"><div></div><div></div><div></div><div></div><div></div><div></div><div></div><div></div></div>
        <div>Processing...</div>
    </div>
    <div id="result-view" hidden>
        <h2>Final output</h2>
        <div id="data-result">
        
        </div>
    </div>
  </div>
`

const startToPeptonize = async function() {
    document.getElementById("result-view").hidden = true;
    document.getElementById("peptonize-button").hidden = true;
    document.getElementById("loading-spinner").hidden = false;
    const data = await (await fetch("data/rescored.psms.tsv")).text();
    document.getElementById("data-result").textContent = await peptonize(data);
    document.getElementById("result-view").hidden = false;
    document.getElementById("loading-spinner").hidden = true;
    // document.getElementById("peptonize-button").hidden = false;
}

document.getElementById("peptonize-button").addEventListener("click", startToPeptonize);
