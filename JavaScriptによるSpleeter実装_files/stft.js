"use strict"

function initWindow(frame_size) {
  const win = [];
  for (var i = 0; i < frame_size; i++) {
    win[i] = 0.5 - 0.5 * Math.cos(2.0 * Math.PI * i / frame_size);
  }
  return win;
}

// http://yukara-13.hatenablog.com/entry/2013/11/17/210204
function mySTFT(src, frame_size, hop_size, proc_frames) {
  const win = initWindow(frame_size);  // hann窓を作る
  const out = new Array();
  const signal = new Float32Array(frame_size*2);
  
  for (var i = 0; i < proc_frames; i++) {
    const start = hop_size * i - frame_size / 2;
    if (start > src.length) {
      out.push(new Float32Array((frame_size / 2 + 1) * 2));
      continue;
    } else if (start < 0) {
      for (var j = 0; j < frame_size; j++) {
        signal[j*2+0] = src[Math.abs(start + j)] * win[j]; // reflect
        signal[j*2+1] = 0;
      }
    } else if (start + frame_size > src.length) {
      for (var j = 0; j < frame_size; j++) {
        const exc = (start + j) - src.length;
        signal[j*2+0] = ((start + j < src.length) ? src[start + j] : src[src.length - 1 - exc]) * win[j]; // reflect
        signal[j*2+1] = 0;
      }
    } else {
      for (var j = 0; j < frame_size; j++) {
        signal[j*2+0] = src[start + j] * win[j];
        signal[j*2+1] = 0;
      }
    }
    
    const phasors = myFFT(signal);
    const half = new Float32Array((frame_size / 2 + 1)*2);
    for (var j = 0; j < frame_size / 2 + 1; j++) {
      half[j*2+0] = phasors[j*2+0];
      half[j*2+1] = phasors[j*2+1];
    }
    out.push(half);
  }
  return out;
}

function myISTFT(dst, mag, phase, frame_size, hop_size, proc_frames) {
  const win = initWindow(frame_size);  // hann窓を作る
  const fphasors = new Float32Array(frame_size*2);
  for (var i = 0; i < proc_frames; i++) {
    const start = hop_size * i - frame_size / 2;
    if (start > dst.length) {
      return;
    }
    const fmag   = mag  [i];
    const fphase = phase[i];
    
    for (var j = 0; j < frame_size / 2 + 1; j++) {
      fphasors[j*2+0] = Math.cos(fphase[j]) * fmag[j];
      fphasors[j*2+1] = Math.sin(fphase[j]) * fmag[j];
    }
    for (var j = 0; j < frame_size / 2 - 1; j++) {
      fphasors[(frame_size - 1 - j)*2+0] =  fphasors[(1 + j)*2+0];
      fphasors[(frame_size - 1 - j)*2+1] = -fphasors[(1 + j)*2+1];
    }
    
    const signal = myIFFT(fphasors);
    
    if (start < 0) {
      for (var j = -start; j < frame_size; j++) {
        dst[start + j] += signal[j*2+0] * win[j];
      }
    } else if (start + frame_size > dst.length) {
      for (var j = 0; j < frame_size && start + j < dst.length; j++) {
        dst[start + j] += signal[j*2+0] * win[j];
      }
    } else {
      for (var j = 0; j < frame_size; j++) {
        dst[start + j] += signal[j*2+0] * win[j];
      }
    }
  }
}

// http://www.kurims.kyoto-u.ac.jp/~ooura/fftman/ftmn1_2.html
function myFFT(a){
  const n = a.length / 2;
  const theta = -2 * Math.PI / n;
  /* ---- scrambler ---- */
  var i = 0;
  for (var j = 1; j < n - 1; j++) {
    for (var k = n >> 1; k > (i ^= k); k >>= 1);
    if (j < i) {
      var xr = a[j*2+0];
      var xi = a[j*2+1];
      a[j*2+0] = a[i*2+0];
      a[j*2+1] = a[i*2+1];
      a[i*2+0] = xr;
      a[i*2+1] = xi;
    }
  }
  var m;
  for (var mh = 1; (m = mh << 1) <= n; mh = m) {
    var irev = 0;
    for (i = 0; i < n; i += m) {
      var wr = Math.cos(theta * irev);
      var wi = Math.sin(theta * irev);
      for (var k = n >> 2; k > (irev ^= k); k >>= 1);
      for (var j = i; j < mh + i; j++) {
        k = j + mh;
        xr = a[j*2+0] - a[k*2+0];
        xi = a[j*2+1] - a[k*2+1];
        a[j*2+0] += a[k*2+0];
        a[j*2+1] += a[k*2+1];
        a[k*2+0] = wr * xr - wi * xi;
        a[k*2+1] = wr * xi + wi * xr;
      }
    }
  }
  return a;
}

function myIFFT(a) {
  const ca = new Float32Array(a.length);
  for (var i = 0; i < a.length / 2; i++) {
    ca[i*2+0] =  a[i*2+0];
    ca[i*2+1] = -a[i*2+1];
  }
  const ps = myFFT(ca);
  const out = new Float32Array(ps.length);
  for (var i = 0; i < ps.length / 2; i++) {
    out[i*2+0] =  ps[i*2+0] / (ps.length / 2);
    out[i*2+1] = -ps[i*2+1] / (ps.length / 2);
  }
  return out;
}

