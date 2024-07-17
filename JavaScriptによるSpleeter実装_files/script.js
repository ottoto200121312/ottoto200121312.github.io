const MODEL_URL_2 = 'https://d3p65luafbyyca.cloudfront.net/spleeter_2stems_191117/model.json';
const MODEL_URL_4 = 'https://d3p65luafbyyca.cloudfront.net/spleeter_4stems_191117/model.json';
const MODEL_URL_5 = 'https://d3p65luafbyyca.cloudfront.net/spleeter_5stems_191118/model.json';
const FFT_SIZE = 4096;
const PATCH_LENGTH = 512;
const HOP_SIZE = 1024;
const SR = 44100;

var aud = {};

let model;

tf.ENV.set('WEBGL_CONV_IM2COL', false);

function info(msg) {
  document.getElementById('info').innerHTML = msg;
}

function error(msg) {
  document.getElementById('info').innerHTML = '<span class="error-message">' + msg + '</span>';
}

function convert() {
  initialize(false);
}

function initialize() {
  if (aud.dst) {
    const ok = confirm("変換済みデータがありますが、よろしいですか？");
    if (!ok) { return; }
  }
  
  const file = document.getElementById('audioFile').files[0];
  if (!file) {
    info("ファイルを選択してください。"); 
    setControlsDefault();
    return;
  }
  
  setControlsProcessing();
  
  const prev_model_url = aud.model_url;
  aud = {};     // GCのために初期化
  aud.model_url = prev_model_url;
  
  // パラメータ読み込み
  const radTypeElems = document.getElementsByName("modelType");
  var i = 0;
  for (; !radTypeElems[i].checked && i < radTypeElems.length; i++) {}
  const modelType = i < radTypeElems.length ? parseInt(radTypeElems[i].value) : 0;
  
  let model_url;
  switch (modelType) {
    case 4:
      model_url = MODEL_URL_4;
      aud.stems = 4;
      break;
    case 5:
      model_url = MODEL_URL_5;
      aud.stems = 5;
      break;
    case 2:       // fall through
    default:
      model_url = MODEL_URL_2;
      aud.stems = 2;
      break;
  }
  
  console.log("url: " + model_url + ", backend: " + tf.getBackend());
  if (aud.model_url == model_url) { // loaded
    info("フレームワークを初期化しています。最大5分程度かかります。操作をせずにお待ちください。");
    setTimeout(function() { warmUp(); }, 10);
  } else{
    info("ニューラルネットワークモデルを読み込んでいます。最大10分程度かかります。操作をせずにお待ちください。");
    model = null;
    tf.loadLayersModel(model_url)
      .then(pretrainedModel => {
        console.log("Model loaded: " + model_url + ", backend: " + tf.getBackend());
        aud.model_url = model_url;
        model = pretrainedModel;
        
        let msg = "フレームワークを初期化しています。最大5分程度かかります。操作をせずにお待ちください。";
        if (tf.getBackend() == "cpu") {
          msg += "　※GPUがサポートされないためCPUを使用しています。処理が極めて遅くなります。";
        }
        info(msg);
        setTimeout(function() { warmUp(); }, 10); // 0だとEdgeで強制リロード
      });
   }
}

function warmUp() {
  const shape = [1, PATCH_LENGTH, FFT_SIZE / 4, 2];
  const x = new Float32Array(1 * PATCH_LENGTH * (FFT_SIZE / 4) * 2); // all zero
  try {
    var tx = tf.tensor(x, shape);
    var ty = predict(tx);
    for (var i = 0; i < aud.stems; i++) { ty[i].dataSync(); }
  } catch (e) {
    console.error(e);
    error("初期化中にエラーが発生しました。恐らくGPUメモリが不足しています。他のアプリ、タブをすべて閉じて再実行するか、他のブラウザ、PCをお試しください。");
    setControlsDefault();
    return;
  } finally {
    tf.dispose(tx);
    for (var i = 0; i < aud.stems; i++) { tf.dispose(ty[i]); }
  }
  info("曲を読み込んでいます。");
  readFile();
}

function predict(input) {
  return tf.tidy(() => {
    return model.predict(input);
  });
}

function setControlsDefault() {
  setConvertControlsEnabled(true);
  setAudioControlsEnabled(false);
}

function setControlsProcessing() {
  setConvertControlsEnabled(false);
  setAudioControlsEnabled(false);
  setStemResults([]);
}

function setControlsCompleted() {
  setConvertControlsEnabled(true);
  setAudioControlsEnabled(true);
}

function setConvertControlsEnabled(en) {
  document.getElementById("btnConvert"     ).disabled = !en;
}

function setAudioControlsEnabled(en) {
  document.getElementById("resultProcessing").className = "row" + (!en ? "" : " row-hide");
  document.getElementById("resultCompleted1").className  = "row" + ( en ? "" : " row-hide");
  document.getElementById("resultCompleted2").className  = "row" + ( en ? "" : " row-hide");
}

function getStemLabels(num_stems) {
  const stem_labels = [['vocals', 'accompaniment'], ["vocals", "drums", "bass", "other"], ["vocals", "piano", "drums", "bass", "other"]];
  return num_stems == 0 ? [] : (stem_labels[(num_stems == 2) ? 0 : (num_stems == 4 ? 1 : 2)]);
}

// [[url, download(label)] , ... ]
function setStemResults(stem_results) {
  const num_stems = stem_results.length;
  const labels = getStemLabels(num_stems);
  
  for (var i = 0; i < 5; i++) {
    const en = (i < num_stems);
    const stem_str = "Stem" + (i+1);
    document.getElementById("para"    + stem_str).className = "four columns " + (en ? "" : "row-hide");
    document.getElementById("spLabel" + stem_str).innerHTML = en ? labels[i] : "";
    
    const audio = document.getElementById("audio"   + stem_str);
    if (en) {
      audio.src = stem_results[i][0];
      audio.id = `audioStem${i+1}`; // IDを設定
    } else {
      audio.removeAttribute('src');
    }
    
    const link = document.getElementById("linkDownload" + stem_str);
    link.href     = en ? stem_results[i][0] : "javascript:void(0)";
    link.download = en ? stem_results[i][1] : null;
  }
}


function readFile() {
  const reader = new FileReader();
  reader.onerror = function() {
    error('読み込み時にエラーが発生しました。コード: ' + reader.error.code);
    setControlsDefault();
  };
 
  reader.onload = function() {
    const arrayBuffer = reader.result;  // Get ArrayBuffer
    info(file.name + "を読み込みました。デコードしています。");
    decode(file.name, arrayBuffer);
  };

  const file = document.getElementById('audioFile').files[0];
  reader.readAsArrayBuffer(file);
}

function decode(fileName, arrayBuffer) {
  const audioCtx = new AudioContext({"sampleRate":SR});
  audioCtx.decodeAudioData(arrayBuffer, 
    function(decodedData) {
      const source = audioCtx.createBufferSource();
      source.buffer = decodedData;
      if (source.buffer.numberOfChannels != 2) {
        info("ステレオのデータのみ対応しています。");
        setControlsDefault();
        return;
      }
      
      resample(fileName, source);
    },
    function(e) {
      error("デコード時にエラーが発生しました。ブラウザがファイルの形式に対応していない可能性があります。他の形式に変換してお試しください。:" + e.name + " " + e.message);
      setControlsDefault();
    });
}

function resample(fileName, bufSrc) {
  const sourceAudioBuffer = bufSrc.buffer;
  // console.log("original SR: " + sourceAudioBuffer.sampleRate);
  if (sourceAudioBuffer.sampleRate == SR) {
    // console.log("buflen: " + sourceAudioBuffer.length + ", dur*sr:" + sourceAudioBuffer.duration * SR);
    info("前処理をしています。");
    aud.src = [sourceAudioBuffer.getChannelData(0), sourceAudioBuffer.getChannelData(1)];
    aud.size = sourceAudioBuffer.length;
    setTimeout(function() { prepare(fileName); }, 10);
    return;
  }
  
  info("リサンプリングしています。");
  const offlineCtx = new OfflineAudioContext(sourceAudioBuffer.numberOfChannels, sourceAudioBuffer.duration * sourceAudioBuffer.numberOfChannels * SR, SR);
  const buffer = offlineCtx.createBuffer(sourceAudioBuffer.numberOfChannels, sourceAudioBuffer.length, sourceAudioBuffer.sampleRate);
  
  // Copy the source data into the offline AudioBuffer
  for (var channel = 0; channel < sourceAudioBuffer.numberOfChannels; channel++) {
    buffer.copyToChannel(sourceAudioBuffer.getChannelData(channel), channel);
  }
  
  // Play it from the beginning.
  const source = offlineCtx.createBufferSource();
  source.buffer = sourceAudioBuffer;
  source.connect(offlineCtx.destination);
  source.start(0);
  
  const validResampledLength = Math.floor(sourceAudioBuffer.duration * SR);
  offlineCtx.oncomplete = function(e) {
    // `resampled` contains an AudioBuffer resampled at SR.
    // use resampled.getChannelData(x) to get an Float32Array for channel x.
    // なぜかreasampledのdurationは倍に伸びている。謎。
    const resampled = e.renderedBuffer;
    const chData0 = resampled.getChannelData(0).slice(0, validResampledLength);
    const chData1 = resampled.getChannelData(1).slice(0, validResampledLength);
    info("前処理をしています。");
    aud.src = [chData0, chData1];
    aud.size = validResampledLength;
    prepare(fileName);
  }
  offlineCtx.startRendering();
}

function prepare(fileName) {
  // 変換後再生用
  const audioBlob  = createWave(aud.src);
  const myURL = window.URL || window.webkitURL;
  const audioR = document.getElementById("audioOriginal");
  const url = myURL.createObjectURL(audioBlob);
  audioR.src = url;
  
  // 処理用のオブジェクトを作る
  aud.numOfPatches = Math.floor(Math.floor((aud.size - 1) / HOP_SIZE) / PATCH_LENGTH) + 1;
  
  // STFTする
  aud.mag   = []; // ch,time,freq
  aud.phase = [];
  aud.mstem = [];
  for (var si = 0; si < aud.stems; si++) { aud.mstem.push([[],[]]); }
  for (var ch = 0; ch < 2; ch++) {
    const spec = mySTFT(aud.src[ch], FFT_SIZE, HOP_SIZE, aud.numOfPatches * PATCH_LENGTH);
    
    const mag   = [];
    const phase = [];
    var mmax = 0;
    for (var i = 0; i < spec.length; i++) {
      const fmag   = new Float32Array(FFT_SIZE / 2 + 1);
      const fphase = new Float32Array(FFT_SIZE / 2 + 1);
      for (var j = 0; j < FFT_SIZE / 2 + 1; j++) {
        fmag  [j] = Math.sqrt(Math.pow(spec[i][j*2+1], 2) + Math.pow(spec[i][j*2+0], 2));
        fphase[j] = Math.atan2(spec[i][j*2+1], spec[i][j*2+0]);
      }
      mag  .push(fmag);
      phase.push(fphase);
    }
    
    aud.mag.push(mag);
    aud.phase.push(phase);
  }
  aud.src = null; // GC
  
  const pos = fileName.lastIndexOf('.');
  aud.base_name = fileName.slice(0, pos);
  
  info("処理を開始しています。");
  setTimeout(function() { processPatch(0); }, 10);
}

function processPatch(pIndex) {
  if (pIndex * PATCH_LENGTH * HOP_SIZE + FFT_SIZE > aud.size) {
    info("後処理をしています。");
    setTimeout(function() { postProcess(); }, 10);
    return;
  }
  
  const success = inference(aud.mag, aud.mstem, pIndex * PATCH_LENGTH);
  if (!success) {
    error("処理中にエラーが発生しました。恐らくGPUメモリが不足しています。他のアプリ、タブをすべて閉じて再実行するか、他のブラウザ、PCをお試しください。");
    setControlsDefault();
    return;
  }
  
  info((pIndex + 1) + "/" + aud.numOfPatches + "を処理しています。");
  setTimeout(function() { processPatch(pIndex + 1); }, 10);
}

function inference(mag, mstem, magIndex) {
  const EPS = 1e-10;
  const INF_FREQ = FFT_SIZE / 4;
  const PATCH_SIZE = 1 * PATCH_LENGTH * INF_FREQ * 2;
  
  // magを配列に設定する
  const x = new Float32Array(PATCH_SIZE); // [time,freq,ch]
  let inpMagAllZero = true;
  for (var i = 0; i < INF_FREQ; i++) {
    for (var j = 0; j < PATCH_LENGTH; j++) {
      const xi = (j * INF_FREQ + i) * 2;
      x[xi + 0] = mag[0][magIndex + j][i];
      x[xi + 1] = mag[1][magIndex + j][i];
      inpMagAllZero &= (x[xi + 0] == 0 && x[xi + 1] == 0);
    }
  }
  
  // inferenceする
  let outMagAllZero = true;
  const stems = aud.stems;
  const shape = [1, PATCH_LENGTH, INF_FREQ, 2];
  try {
    var tx = tf.tensor(x, shape);
    var ty = predict(tx);
    
    const y = [];
    for (var si = 0; si < stems; si++) { y.push(ty[si].dataSync()); }
    
    for (var i = 0; i < PATCH_SIZE; i++) {
      let sum = EPS;
      for (var si = 0; si < stems; si++) {
        const v = y[si][i];
        sum += v * v;
      }
      for (var si = 0; si < stems; si++) {
        const v = y[si][i];
        y[si][i] = (v * v + (EPS / stems)) / sum;
        outMagAllZero &= (v == 0);
      }
    }
    
    for (var si = 0; si < stems; si++) {
      for (var i = 0; i < PATCH_LENGTH; i++) {
        fsmag0 = new Float32Array(INF_FREQ);
        fsmag1 = new Float32Array(INF_FREQ);
        for (var j = 0; j < INF_FREQ; j++) {
          const yi = (i * INF_FREQ + j) * 2;
          const y0 = y[si][yi + 0]
          const y1 = y[si][yi + 1];
          fsmag0[j] = y0;
          fsmag1[j] = y1;
        }
        mstem[si][0].push(fsmag0);
        mstem[si][1].push(fsmag1);
      }
    }
  } catch (error) {
    console.error(error);
    return false;
  } finally {
    try {
      tf.dispose(tx);
      for (var i = 0; i < aud.stems; i++) { tf.dispose(ty[i]); }
    } catch (err2) {
      console.error(err2);
    }
  }
  if (outMagAllZero) { console.log("Inference result is null: " + magIndex); }
  return inpMagAllZero || !outMagAllZero; // true for success
}

function postProcessMag(mstem, morg, phase) {      // magを計算し、帯域を拡張する
  const nmag   = [[],[]]
  const nphase = [[],[]]
  for (var ch = 0; ch < 2; ch++) {
    for (var i = 0; i < morg[0].length; i++) {
      const fnmag   = new Float32Array(FFT_SIZE / 2 + 1);
      const fnphase = new Float32Array(FFT_SIZE / 2 + 1);
      var j = 0;
      for (; j < FFT_SIZE / 4; j++) {
        fnmag  [j] = mstem[ch][i][j] * morg[ch][i][j];
        fnphase[j] = phase[ch][i][j];
      }
      for (; j < FFT_SIZE / 2 + 1; j++) {
        fnmag  [j] = 0;
        fnphase[j] = phase[ch][i][j];
      }
      nmag  [ch].push(fnmag);
      nphase[ch].push(fnphase);
    }
  }
  return [nmag, nphase];
}

function downSampleWaveform(waveform, targetLength) {
  let result = [];
  let stepSize = Math.floor(waveform.length / targetLength);
  for (let i = 0; i < targetLength; i++) {
    result.push(waveform[i * stepSize]);
  }
  return result;
}

function postProcess() {
  const stems = aud.stems;
  
  aud.dst = []
  for (var si = 0; si < stems; si++) {
    aud.dst.push([new Float32Array(aud.size), new Float32Array(aud.size)]);
  }
  
  // ISTFTする
  for (var si = 0; si < stems; si++) {
    const mp = postProcessMag(aud.mstem[si], aud.mag, aud.phase);
    myISTFT(aud.dst[si][0], mp[0][0], mp[1][0], FFT_SIZE, HOP_SIZE, aud.numOfPatches * PATCH_LENGTH);
    myISTFT(aud.dst[si][1], mp[0][1], mp[1][1], FFT_SIZE, HOP_SIZE, aud.numOfPatches * PATCH_LENGTH);
  }
  aud.mstem = null;     // gc
  
  // clipする : hann窓+HOPで1.5倍になっている
  for (var si = 0; si < stems; si++) {
    for (var ch = 0; ch < 2; ch++) {
      const dst = aud.dst[si][ch];
      for (var i = 0; i < aud.size; i++) {
        dst[i] /= 1.5;
        if      (dst[i] >  1) { dst[i] =  1; }
        else if (dst[i] < -1) { dst[i] = -1; }
      }
    }
  }
  
  const results = []
  const labels  = getStemLabels(aud.stems);
  for (var si = 0; si < stems; si++) {
    const audioBlob  = createWave(aud.dst[si]);
    const myURL = window.URL || window.webkitURL;
    const url  = myURL.createObjectURL(audioBlob);
    results.push([url, aud.base_name + "_" + labels[si] + ".wav"]);
  }
  setStemResults(results);
  
  setControlsCompleted();
  info("変換が完了しました。");

  waveforms = [];
  for (var si = 0; si < stems; si++) {
    let waveform = [];
    for (let i = 0; i < aud.size; i++) {
      waveform.push((aud.dst[si][0][i] + aud.dst[si][1][i]) / 2);
    }
    let downSampledWaveform = downSampleWaveform(waveform, 1000);
    waveforms.push(downSampledWaveform);
  }

 // p5.jsのセットアップと描画を開始
 new p5(function(p) {
  p.setup = function() {
    setupP5(p);
  };
  p.draw = function() {
    drawP5(p);
  };
}, 'waveformCanvas');
}



// export to wav
function createWave(dst) {
  const length = dst[0].length * 2 * 2 + 44;
  const buffer = new ArrayBuffer(length);
  const view = new DataView(buffer);
  let pos = 0;

  // write WAVE header
  setUint32(0x46464952);                         // "RIFF"
  setUint32(length - 8);                         // file length - 8
  setUint32(0x45564157);                         // "WAVE"

  setUint32(0x20746d66);                         // "fmt " chunk
  setUint32(16);                                 // length = 16
  setUint16(1);                                  // PCM (uncompressed)
  setUint16(2); // numOfChan);
  setUint32(SR); // abuffer.sampleRate);
  setUint32(SR * 2 * 2); //abuffer.sampleRate * 2 * numOfChan); // avg. bytes/sec
  setUint16(2 * 2); // numOfChan * 2);                      // block-align
  setUint16(16);                                 // 16-bit (hardcoded in this demo)

  setUint32(0x61746164);                         // "data" - chunk
  setUint32(length - pos - 4);                   // chunk length

  // write interleaved data
  let offset = 0;
  while(pos < length) {
    for(ch = 0; ch < 2; ch++) {             // interleave channels
      sample = Math.max(-1, Math.min(1, dst[ch][offset])); // 必要ないはずだけど念のため clamp
      sample = ((sample < 0) ? sample * 0x8000 : sample * 0x7fff) | 0; // scale to 16-bit signed int
      view.setInt16(pos, sample, true);          // write 16-bit sample
      pos += 2;
    }
    offset++;                                     // next source sample
  }
  // create Blob
  return new Blob([buffer], {type: "audio/wav"});

  function setUint16(data) {
    view.setUint16(pos, data, true);
    pos += 2;
  }

  function setUint32(data) {
    view.setUint32(pos, data, true);
    pos += 4;
  }
}

document.getElementById("btnConvert"     ).addEventListener("click", convert);

setControlsDefault();

window.onerror = function(msg, url, line, col, error) {
  //  https://stackoverflow.com/questions/951791/javascript-global-error-handling
  console.log(error);
  var extra = !col ? '' : ', column: ' + col;
  extra += !error ? '' : ', error: ' + error;
  
  error("エラーが発生しました。ブラウザを閉じて開きなおすと改善する可能性があります。\n" + msg +" (" + line + extra + ")");
  setControlsDefault();
  
  var suppressErrorAlert = true;
  return suppressErrorAlert;
};



// <<p5
let waveforms = [];
let audioContexts = [];
let analyserNodes = [];
let sourceNodes = [];
let isPlaying = false;
let playAllButton;
let timelineSlider;
let currentTimeSpan;
let totalTimeSpan;
let duration = 0;
let drumPeaks = [];
let drumLagDuration = 5; // フレーム数でラグの長さを指定
let pianoNotes = [];
let pianoDelay = 2; // ピアノのディレイフレーム数



function setupP5(p) {
  let canvas = p.createCanvas(800, 400);
  canvas.parent('waveformCanvas');

  // 各stemのオーディオコンテキストとアナライザーノードを設定
  for (let i = 0; i < aud.stems; i++) {
    let audio = document.getElementById(`audioStem${i+1}`);
    if (audio) {
      let audioContext = new (window.AudioContext || window.webkitAudioContext)();
      let analyser = audioContext.createAnalyser();
      analyser.fftSize = 2048;
      
      let source = audioContext.createMediaElementSource(audio);
      source.connect(analyser);
      analyser.connect(audioContext.destination);

      audioContexts.push(audioContext);
      analyserNodes.push(analyser);
      sourceNodes.push(source);

      // 最長のオーディオの長さを取得
      duration = Math.max(duration, audio.duration);
    }
  }

  // コントロール要素の設定
  playAllButton = document.getElementById('playAllButton');
  playAllButton.addEventListener('click', togglePlayAll);

  timelineSlider = document.getElementById('timelineSlider');
  timelineSlider.max = duration;
  timelineSlider.addEventListener('input', seekAudio);

  currentTimeSpan = document.getElementById('currentTime');
  totalTimeSpan = document.getElementById('totalTime');
  totalTimeSpan.textContent = formatTime(duration);

  // タイマーを設定して現在の再生時間を更新
  setInterval(updateCurrentTime, 1000);
}

function togglePlayAll() {
  let audios = [];
  for (let i = 0; i < aud.stems; i++) {
    let audio = document.getElementById(`audioStem${i+1}`);
    if (audio) {
      audios.push(audio);
    }
  }

  if (!isPlaying) {
    // すべてのオーディオを同時に再生開始
    audios.forEach(audio => audio.play());
    isPlaying = true;
    playAllButton.textContent = 'Pause All';
  } else {
    // すべてのオーディオを停止
    audios.forEach(audio => audio.pause());
    isPlaying = false;
    playAllButton.textContent = 'Play All';
  }
}

function seekAudio() {
  let time = timelineSlider.value;
  for (let i = 0; i < aud.stems; i++) {
    let audio = document.getElementById(`audioStem${i+1}`);
    if (audio) {
      audio.currentTime = time;
    }
  }
  updateCurrentTime();
}

function updateCurrentTime() {
  let currentTime = 0;
  if (aud.stems > 0) {
    let firstAudio = document.getElementById('audioStem1');
    if (firstAudio) {
      currentTime = firstAudio.currentTime;
    }
  }
  currentTimeSpan.textContent = formatTime(currentTime);
  timelineSlider.value = currentTime;
}

function formatTime(time) {
  let minutes = Math.floor(time / 60);
  let seconds = Math.floor(time % 60);
  return `${minutes}:${seconds.toString().padStart(2, '0')}`;
}

function drawP5(p) {
  p.background(240);
  
  const stemLabels = getStemLabels(aud.stems);
  
  for (let i = 0; i < analyserNodes.length; i++) {
    let y = i * p.height / analyserNodes.length;
    let h = p.height / analyserNodes.length;
    
    // 波形データを取得
    let bufferLength = analyserNodes[i].frequencyBinCount;
    let dataArray = new Uint8Array(bufferLength);
    analyserNodes[i].getByteTimeDomainData(dataArray);
    
    p.stroke(0);
    p.noFill();
    
    if (stemLabels[i].toLowerCase() === 'drums') {
      // ドラムの処理（前回のコードと同じ）
      p.beginShape();
      let isPeak = false;
      let peakValue = 0;
      let peakIndex = 0;
      
      for (let j = 0; j < bufferLength; j++) {
        let x = p.map(j, 0, bufferLength, 0, p.width);
        let v = p.map(dataArray[j], 0, 255, 0, h);
        
        if (v > h/2 + h/4 && v > peakValue) {
          isPeak = true;
          peakValue = v;
          peakIndex = j;
        }
        
        if (isPeak && j > peakIndex + 5) {
          isPeak = false;
          drumPeaks.push({x: x, y: y + peakValue, frame: p.frameCount});
        }
        
        p.vertex(x, y + v);
      }
      p.endShape();
      
      p.fill(255, 0, 0);
      p.noStroke();
      for (let k = drumPeaks.length - 1; k >= 0; k--) {
        let peak = drumPeaks[k];
        let age = p.frameCount - peak.frame;
        if (age < drumLagDuration) {
          p.ellipse(peak.x, peak.y, 10, 10);
        } else {
          drumPeaks.splice(k, 1);
        }
      }
    } else if (stemLabels[i].toLowerCase() === 'bass') {
      // ベースの処理（前回のコードと同じ）
      p.beginShape(p.LINES);
      for (let j = 0; j < bufferLength; j += 2) {
        let x1 = p.map(j, 0, bufferLength, 0, p.width);
        let x2 = p.map(j + 1, 0, bufferLength, 0, p.width);
        let v = p.map(dataArray[j], 0, 255, 0, h);
        p.vertex(x1, y + h/2);
        p.vertex(x1, y + v);
        p.vertex(x1, y + v);
        p.vertex(x2, y + v);
        p.vertex(x2, y + v);
        p.vertex(x2, y + h/2);
      }
      p.endShape();
    } else if (stemLabels[i].toLowerCase() === 'piano') {
      // ピアノの処理
      drawPianoStaff(p, y, h, dataArray, bufferLength);
    } else {
      // その他のステムは通常の波形を描画
      p.beginShape();
      for (let j = 0; j < bufferLength; j++) {
        let x = p.map(j, 0, bufferLength, 0, p.width);
        let v = p.map(dataArray[j], 0, 255, 0, h);
        p.vertex(x, y + v);
      }
      p.endShape();
    }
    
    // ステムラベルを表示
    p.fill(0);
    p.noStroke();
    p.textAlign(p.LEFT, p.TOP);
    p.text(stemLabels[i], 10, y + 10);
  }
  // タイムラインマーカーを描画
  let markerX = p.map(timelineSlider.value, 0, duration, 0, p.width);
  p.stroke(255, 0, 0);
  p.line(markerX, 0, markerX, p.height);
}

function drawPianoStaff(p, y, h, dataArray, bufferLength) {
  // 五線譜を描画
  p.stroke(0);
  for (let i = 0; i < 5; i++) {
    let lineY = y + h * (0.3 + i * 0.1);
    p.line(0, lineY, p.width, lineY);
  }
  
  // ピアノの音符を描画
  for (let j = 0; j < bufferLength; j += 10) {
    let x = p.map(j, 0, bufferLength, 0, p.width);
    let v = p.map(dataArray[j], 0, 255, 0, h);
    if (v > h/2) {
      let noteY = y + h * (1 - v/h);
      pianoNotes.push({x: x, y: noteY, frame: p.frameCount});
    }
  }
  
  // 音符とディレイを描画
  p.noStroke();
  for (let k = pianoNotes.length - 1; k >= 0; k--) {
    let note = pianoNotes[k];
    let age = p.frameCount - note.frame;
    if (age < pianoDelay) {
      // メイン音符（四角形）
      p.fill(0);
      p.rect(note.x, note.y - 4, 6, 8);
      
      // ディレイエフェクト（薄くなる四角形）
      // for (let d = 1; d <= 2; d++) {
      //   let delayAlpha = p.map(d, 1, 2, 150, 50);
      //   p.fill(0, delayAlpha);
      //   p.rect(note.x - d * 8, note.y - 4, 6, 8);
      // }
    } else {
      pianoNotes.splice(k, 1);
    }
  }
}




