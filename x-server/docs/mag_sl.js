// MAG_form fields enabling & disabling + fields validation
// Version-2.0 2013/05/19  fixed setOnloads as it did not work properly in IE after hitting "back" button

function LargeFloat(webpage, windowname, pagewidth, pageheight) {
  launch = window.open (webpage, windowname, 'width='+pagewidth+',height='+pageheight
                       +',toolbar=0,location=0,directories=0,status=0,menubar=0,scrollbars=1,resizable=1');
}

function setOnloads() {
    switchXrayFields();
    switchPolarization();
    switchSubstrate();
}

function getSubway() {
  for (var i=0; i < document.magform.subway.length; i++) {
    if (document.magform.subway[i].checked) {
      return document.magform.subway[i].value;
    }
  }
  return -1;
}

function switchXrayFields() {
  if (document.magform.xway.value != 3) {
    document.magform.wave.disabled=false;
    document.magform.line.disabled=true;
  } else {
    document.magform.wave.disabled=true;
    document.magform.line.disabled=false;
  }
}

function switchPolarization() {
  if (document.magform.ipol.value != 3) {
    document.magform.polangle.disabled=true;
  } else {
    document.magform.polangle.disabled=false;
  }
}

function switchSubstrate() {
  var subway = getSubway();
  if (subway == 1) {
    document.magform.code.disabled=false;
    document.magform.chem.disabled=true;
    document.magform.rho.disabled=true;
    document.magform.x0.disabled=true;
  }
  else if (subway == 2) {
    document.magform.code.disabled=true;
    document.magform.chem.disabled=false;
    document.magform.rho.disabled=false;
    document.magform.x0.disabled=true;
  }
  else if (subway == 3) {
    document.magform.code.disabled=true;
    document.magform.chem.disabled=true;
    document.magform.rho.disabled=true;
    document.magform.x0.disabled=false;
  }
}

function valid_chars(str) {
// This is full ASCII pattern:
// var pattern = new RegExp(/^[\r\n\s\w!"#$%&'()*+,-./:;<=>?@\\\[\]^{|}~]*$/);
// where \s stands for both spaces and tabs and \w is the same as [a-zA-Z0-9_]
// This is restricted to what we allow in forms (exclude: "$'<>\|)
   var pattern = new RegExp(/^[\r\n\s\w!#%&()*+,-./:;=?@\[\]^{}~]*$/);
   return pattern.test(str);
}

function valid_code_chars(str) {
// This is restricted to what we allow in the material code:
// \w is a metacharacter matching word characters: [a-zA-Z0-9_]
   var pattern = new RegExp(/^[\s\w()+-.]*$/);
   return pattern.test(str);
}

function mag_validate() {
  var wave_to_energy = 12.3981;
  var wave_min = 0.01;
  var wave_max = 1000.;
  var energy_min = (wave_to_energy/wave_max).toFixed(6);
  var energy_max = (wave_to_energy/wave_min).toFixed(1);
  if (document.magform.xway.value == 1) {
    if (document.magform.wave.value == '') {
      alert('X-ray wavelength in Angstrem is not specified');
      return false;
    }
    var wave = parseFloat(document.magform.wave.value);
    if (wave < wave_min || wave > wave_max) {
      alert('X-ray wavelength='+wave+' Angstrom must be in range '+wave_min+'-'+wave_max+' Angstrom');
      return false;
    }
  }
  else if (document.magform.xway.value == 2) {
    if (document.magform.wave.value == '') {
      alert('X-ray energy in KeV is not specified');
      return false;
    }
    var energy = parseFloat(document.magform.wave.value);
    if (energy < energy_min || energy > energy_max) {
      alert('X-ray energy='+energy+' KeV must be in range: '
           +energy_min+'-'+energy_max+' KeV \n'
           +'(wavelength: '+wave_min+'-'+wave_max+' Angstrom)');
      return false;
    }
  }
  else if (document.magform.xway.value == 3) {
    var selectedLine = document.magform.line.selectedIndex;
    if (document.magform.line.options[selectedLine].text == '') {
      alert('X-ray line selected, but not chosen from the list');
      return false;
    }
  }
  var subway = getSubway();
  if (subway == 1) {
    var selectedCode = document.magform.code.selectedIndex;
    if (document.magform.code.options[selectedCode].text == '') {
      alert('Substrate code is not specified');
      return false;
    }
  }
  else if (subway == 2) {
    if (document.magform.chem.value == '') {
      alert('Substrate chemical formula is not specified');
      return false;
    }
    if (document.magform.rho.value == '') {
      alert('Substrate mamagial density is not specified');
      return false;
    }
    if (document.magform.rho.value <= 0.) {
      alert('Substrate mamagial density must be a positive value');
      return false;
    }
  }
  else if (subway == 3) {
    var x0string = document.magform.x0.value;
    if (x0string == '') {
      alert('Substrate X0 is not specified');
      return false;
    }
    var x0array = x0string.split(/[, ]/);	//expexted to be array of 2
    if (x0array.length != 2) {
      alert('Substrate X0 must contain X0r and X0i while '+x0array.length+' element(s) found: '+x0array);
      return false;
    }
  }
  if (document.magform.w0.value <= 0 || document.magform.w0.value > 99) {
    alert('The Debye-Waller modifier of x0 should be in range [0-99], w0>0');
    return false;
  }
  if (document.magform.sigma.value < 0) {
    alert('Interface roughness sigma cannot be negative');
    return false;
  }
  if (document.magform.tr.value < 0) {
    alert('Transition layer thickness cannot be negative');
    return false;
  }
  if (document.magform.sigma.value > 0 && document.magform.tr.value > 0) {
    alert('Interface roughness and transition layer have the same effect and cannot be used simultaneously');
    return false;
  }
  var factor = 1.;
  if      (document.magform.unis.value == 0) {factor = 1.;}		//degr.
  else if (document.magform.unis.value == 1) {factor = 60.;}		//min.
  else if (document.magform.unis.value == 2) {factor = 17.4532925;}	//mrad
  else if (document.magform.unis.value == 3) {factor = 3600;}		//sec
  else if (document.magform.unis.value == 4) {factor = 17453.2925;}	//urad
  var scanmin = parseFloat(document.magform.scanmin.value) / factor;
  var scanmax = parseFloat(document.magform.scanmax.value) / factor;
  if (isNaN(scanmin)) {
    alert('Minimum scan range is not numeric');
    return false;
  }
  if (isNaN(scanmax)) {
    alert('Maximum scan range is not numeric');
    return false;
  }
  if (document.magform.unis.value < 5) {				//angular scan
    if (scanmin < 0 || scanmin > 90 ||
        scanmin < 0 || scanmin > 90) {
      alert('Incidence angle must be in range [0-90] degr.');
      return false;
    }
  } else {								//qz-scan
    if (scanmin < 0 || scanmax < 0) {
      alert('qz scan parameter cannot be negative');
      return false;
    }
  }
  var nscanMax = 10001;
  if (document.magform.nscan.value < 1 || document.magform.nscan.value > nscanMax) {
    alert('Number of scan points must be in range [1-'+nscanMax+']');
    return false;
  }

  if (! valid_chars(document.magform.profile.value)) {
    alert('Illegal characters in the surface profile input.\n'+
          'Illegal characters are non-English letters, any\n' +
          'characters encoded in UTF-8, and "$\'<>|\\ \n' +
          'Possibly try to disable browser translating.');
    return false;
  }

  if (! valid_code_chars(document.magform.code.value)) {
    alert('Illegal characters in the crystal/material code input:\n\n'+
          document.magform.code.value + '\n\n' +
          'Illegal characters are non-English letters, any\n' +
          'characters encoded in UTF-8, and the special symbols\n' +
          'other than "()+-.". Possibly try to disable translating\n' +
          'into your language in the browser.');
    return false;
  }

// alert('OK '); return false;
   return true;
}
