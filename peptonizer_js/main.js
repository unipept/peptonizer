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
  </div>
`

const data = await (await fetch("data/rescored.psms.tsv")).text();
await peptonize(data);
