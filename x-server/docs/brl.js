// BRL_form fields enabling & disabling + fields validation
// Version-2.0 2013/05/19  fixed setOnloads as it did not work properly in IE after hitting "back" button

function LargeFloat(webpage, windowname, pagewidth, pageheight) {
  launch = window.open (webpage, windowname, 'width='+pagewidth+',height='+pageheight
                       +',toolbar=0,location=0,directories=0,status=0,menubar=0,scrollbars=1,resizable=1');
}

function setOnloads() {
    switchXrayFields();
}

function getXway() {
  for (var i=0; i < document.brlform.xway.length; i++) {
    if (document.brlform.xway[i].checked) {
      return document.brlform.xway[i].value;
    }
  }
  return -1;
}

function switchXrayFields() {
  var xway = getXway();
  if (xway == 1 || xway == 2) {
    document.brlform.wave.disabled=false;
    document.brlform.line.disabled=true;
    document.brlform.i31.disabled=true;
    document.brlform.i32.disabled=true;
    document.brlform.i33.disabled=true;
  }
  else if (xway == 3) {
    document.brlform.wave.disabled=true;
    document.brlform.line.disabled=false;
    document.brlform.i31.disabled=true;
    document.brlform.i32.disabled=true;
    document.brlform.i33.disabled=true;
  }
  else if (xway == 4) {
    document.brlform.wave.disabled=true;
    document.brlform.line.disabled=true;
    document.brlform.i31.disabled=true;
    document.brlform.i32.disabled=true;
    document.brlform.i33.disabled=true;
  }
  else if (xway == 5) {
    document.brlform.wave.disabled=true;
    document.brlform.line.disabled=true;
    document.brlform.i31.disabled=false;
    document.brlform.i32.disabled=false;
    document.brlform.i33.disabled=false;
  }
  else {
    document.brlform.wave.disabled=false;
    document.brlform.line.disabled=false;
    document.brlform.i31.disabled=false;
    document.brlform.i32.disabled=false;
    document.brlform.i33.disabled=false;
  }
}

function brl_validate() {
  var wave_to_energy = 12.3981;
  var wave_min       = 0.1;
  var wave_max       = 10.;
  var energy_min = (wave_to_energy/wave_max).toFixed(4);
  var energy_max = (wave_to_energy/wave_min).toFixed(2);
  var xway = getXway();
  if (xway == 1) {
    var wave = parseFloat(document.brlform.wave.value);
    if (wave < wave_min || wave > wave_max) {
      alert('X-ray wavelength='+wave+' Angstrom must be in range '+wave_min+'-'+wave_max+' Angstrom');
      return false;
    }
  }
  else if (xway == 2) {
    var energy = parseFloat(document.brlform.wave.value);
    if (energy < energy_min || energy > energy_max) {
      alert('X-ray energy='+energy+' KeV must be in range: '
           +energy_min.toFixed(4)+'-'+energy_max.toFixed(2)+' KeV \n'
           +'(wavelength: '+wave_min+'-'+wave_max+' Angstrom)');
      return false;
    }
  }
  else if (xway == 3) {
    var selectedLine = document.brlform.line.selectedIndex;
    if (document.brlform.line.options[selectedLine].text == '') {
      alert('X-ray line selected, but not chosen from the list');
      return false;
    }
  }
  else if (xway == 5) {
    if ((document.brlform.i31.value == 0 || document.brlform.i31.value.length == 0) &&
        (document.brlform.i32.value == 0 || document.brlform.i32.value.length == 0) &&
        (document.brlform.i33.value == 0 || document.brlform.i33.value.length == 0)) {
      alert('Please specify non-zero indices for Bragg reflection-3');
      return false;
    }
  }
  if ((document.brlform.i11.value == 0 || document.brlform.i11.value.length == 0) &&
      (document.brlform.i12.value == 0 || document.brlform.i12.value.length == 0) &&
      (document.brlform.i13.value == 0 || document.brlform.i13.value.length == 0)) {
    alert('Please specify non-zero indices for Bragg reflection-1');
    return false;
  }
  if ((document.brlform.i21.value == 0 || document.brlform.i21.value.length == 0) &&
      (document.brlform.i22.value == 0 || document.brlform.i22.value.length == 0) &&
      (document.brlform.i23.value == 0 || document.brlform.i23.value.length == 0)) {
    alert('Please specify non-zero indices for Bragg reflection-2');
    return false;
  }
// alert('OK ' + document.brlform.nscan.value); return false;
   return true;
}
