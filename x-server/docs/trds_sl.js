// TRDS_form fields enabling & disabling + fields validation
// Version-2.0 2013/05/19  fixed setOnloads as it did not work properly in IE after hitting "back" button

function LargeFloat(webpage, windowname, pagewidth, pageheight) {
  launch = window.open (webpage, windowname, 'width='+pagewidth+',height='+pageheight
                       +',toolbar=0,location=0,directories=0,status=0,menubar=0,scrollbars=1,resizable=1');
}

function setOnloads() {
    switchXrayFields();
    switchSubstrate();
    switchScan();
    switchModel();
}

function getFactor(x) {
  var factor = 1;
  if      (x.value == 0) {factor = 1.;}			//degr.
  else if (x.value == 1) {factor = 60.;}		//min.
  else if (x.value == 2) {factor = 17.4532925;}		//mrad
  else if (x.value == 3) {factor = 3600.;}		//sec
  else if (x.value == 4) {factor = 17453.2925;}		//urad
  return factor;
}

function getSubway() {
  for (var i=0; i < document.trdsform.subway.length; i++) {
    if (document.trdsform.subway[i].checked) {
      return document.trdsform.subway[i].value;
    }
  }
  return -1;
}

function getModel() {
  for (var i=0; i < document.trdsform.model.length; i++) {
    if (document.trdsform.model[i].checked) {
      return document.trdsform.model[i].value;
    }
  }
  return -1;
}

function switchXrayFields() {
  if (document.trdsform.xway.value != 3) {
    document.trdsform.wave.disabled=false;
    document.trdsform.line.disabled=true;
  } else {
    document.trdsform.wave.disabled=true;
    document.trdsform.line.disabled=false;
  }
}

function switchSubstrate() {
  var subway = getSubway();
  if (subway == 1) {
    document.trdsform.code.disabled=false;
    document.trdsform.chem.disabled=true;
    document.trdsform.rho.disabled=true;
    document.trdsform.x0.disabled=true;
  }
  else if (subway == 2) {
    document.trdsform.code.disabled=true;
    document.trdsform.chem.disabled=false;
    document.trdsform.rho.disabled=false;
    document.trdsform.x0.disabled=true;
  }
  else if (subway == 3) {
    document.trdsform.code.disabled=true;
    document.trdsform.chem.disabled=true;
    document.trdsform.rho.disabled=true;
    document.trdsform.x0.disabled=false;
  }
}

function switchScan() {
  if (document.trdsform.scan.value != 4) {
    document.trdsform.unis.disabled=false;
    document.trdsform.unia.disabled=true;
  } else {
    document.trdsform.unis.disabled=true;
    document.trdsform.unia.disabled=false;
  }
}

function switchModel() {
  var model = getModel();
  if (model == 1 || model == 2) {			//uncorrelated or fully correlated
    document.trdsform.Lv.disabled=true;
  } else {
    document.trdsform.Lv.disabled=false;
  }
  if (model == 1) {					//uncorrelated
    document.trdsform.skew.disabled=true;
    document.trdsform.uniw.disabled=true;
  } else {
    document.trdsform.skew.disabled=false;
    document.trdsform.uniw.disabled=false;
  }
  if (model == 4) {					//Lagally
    document.trdsform.Lh2.disabled=false;
  } else {
    document.trdsform.Lh2.disabled=true;
  }
  if (model == 8) {					//Smoothed Pikute
    document.trdsform.stepH.disabled=false;
  } else {
    document.trdsform.stepH.disabled=true;
  }
  if (model == 9) {					//Pershan
    document.trdsform.spread.disabled=false;
  } else {
    document.trdsform.spread.disabled=true;
  }
  if (model < 7) {					//any non-stepped
    document.trdsform.jagg.disabled=false;
    document.trdsform.Qm.disabled=true;
    document.trdsform.unim.disabled=true;
    document.trdsform.add.disabled=true;
  } else {
    document.trdsform.jagg.disabled=true;
    document.trdsform.Qm.disabled=false;
    document.trdsform.unim.disabled=false;
    document.trdsform.add.disabled=false;
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

function trds_validate() {
  var wave_to_energy = 12.3981;
  var wave_min = 0.1;
  var wave_max = 10.;
  var energy_min = (wave_to_energy/wave_max).toFixed(4);
  var energy_max = (wave_to_energy/wave_min).toFixed(2);
  if (document.trdsform.xway.value == 1) {
    if (document.trdsform.wave.value == '') {
      alert('X-ray wavelength in Angstrem is not specified');
      return false;
    }
    var wave = parseFloat(document.trdsform.wave.value);
    if (wave < wave_min || wave > wave_max) {
      alert('X-ray wavelength='+wave+' Angstrom must be in range '+wave_min+'-'+wave_max+' Angstrom');
      return false;
    }
  }
  else if (document.trdsform.xway.value == 2) {
    if (document.trdsform.wave.value == '') {
      alert('X-ray energy in KeV is not specified');
      return false;
    }
    var energy = parseFloat(document.trdsform.wave.value);
    if (energy < energy_min || energy > energy_max) {
      alert('X-ray energy='+energy+' KeV must be in range: '
           +energy_min+'-'+energy_max+' KeV \n'
           +'(wavelength: '+wave_min+'-'+wave_max+' Angstrom)');
      return false;
    }
  }
  else if (document.trdsform.xway.value == 3) {
    var selectedLine = document.trdsform.line.selectedIndex;
    if (document.trdsform.line.options[selectedLine].text == '') {
      alert('X-ray line selected, but not chosen from the list');
      return false;
    }
  }
  var subway = getSubway();
  if (subway == 1) {
    var selectedCode = document.trdsform.code.selectedIndex;
    if (document.trdsform.code.options[selectedCode].text == '') {
      alert('Substrate code is not specified');
      return false;
    }
  }
  else if (subway == 2) {
    if (document.trdsform.chem.value == '') {
      alert('Substrate chemical formula is not specified');
      return false;
    }
    if (document.trdsform.rho.value == '') {
      alert('Substrate matrdsial density is not specified');
      return false;
    }
    if (document.trdsform.rho.value <= 0.) {
      alert('Substrate matrdsial density must be a positive value');
      return false;
    }
  }
  else if (subway == 3) {
    var x0string = document.trdsform.x0.value;
    if (x0string == '') {
      alert('Substrate X0 is not specified');
      return false;
    }
    var x0array = x0string.split(/[, ]/);		//expected to be array of 2
    if (x0array.length != 2) {
      alert('Substrate X0 must contain X0r and X0i while '+x0array.length+' element(s) found: '+x0array);
      return false;
    }
  }
  if (document.trdsform.w0.value <= 0 || document.trdsform.w0.value > 99) {
    alert('The Debye-Waller modifier of x0 should be in range [0-99], w0>0');
    return false;
  }
  if (document.trdsform.sigma.value < 0) {
    alert('Interface roughness sigma cannot be negative');
    return false;
  }
  if (document.trdsform.Lh.value < 0) {
    alert('Lateral correlation length should be positive');
    return false;
  }

  var model = getModel();
  if (model < 7) {
    if (parseFloat(document.trdsform.jagg.value) < 0.1 ||
        parseFloat(document.trdsform.jagg.value) > 1.0) {
      alert('Roughness jaggness should be in range [0.1-1.0]');
      return false;
    }
  }
  if (model == 4) {
    if (parseFloat(document.trdsform.Lh2.value) < parseFloat(document.trdsform.Lh.value)) {
      alert('Lagallys model cross correlation length cannot be less than lateral correlation length');
      return false;
    }
  }
  var angle = 1;
  if (model != 1) {
    angle = parseFloat(document.trdsform.skew.value) / getFactor(document.trdsform.uniw);
    if (angle  > 89. || angle < -89.) {
      alert('Angle of skew roughness transfer should be in range [-89:89] degr.');
      return false;
    }
  }
  if (model >= 7) {
    angle = parseFloat(document.trdsform.Qm.value) / getFactor(document.trdsform.unim);
    if (angle  > 20. || angle < -20.) {
      alert('Stepped surface miscut angle should be in range [-20:20] degr.');
      return false;
    }
  }
  if (model == 8) {
    if (document.trdsform.stepH.value < 0) {
      alert('RMS height of steps in the Pikute model cannot be negative');
      return false;
    }
  }
  if (model == 9) {
    if (parseFloat(document.trdsform.spread.value) < 0 ||
        parseFloat(document.trdsform.spread.value) > parseFloat(document.trdsform.Lh.value)) {
      alert('Terrace spread in the Pershan model should be in range [0:'+document.trdsform.Lh.value+']');
      return false;
    }
  }
  var nscanMax = 1001;
  if (parseFloat(document.trdsform.nscan.value) < 1 ||
      parseFloat(document.trdsform.nscan.value) > nscanMax) {
    alert('Number of scan points must be in range [1-'+nscanMax+']');
    return false;
  }
  if (parseFloat(document.trdsform.noff.value) < 0 ||
      parseFloat(document.trdsform.noff.value) > nscanMax) {
    alert('Number of scan offsets must be in range [0-'+nscanMax+']');
    return false;
  }
  var factor = 1.;
  if (document.trdsform.scan.value != 4) {factor = getFactor(document.trdsform.unis);}	//not a qz,qx scan
  var scanmin = document.trdsform.scanmin.value / factor;
  var scanmax = document.trdsform.scanmax.value / factor;
  var offmin  = document.trdsform.offmin.value  / factor;
  var offmax  = document.trdsform.offmax.value  / factor;
  if (isNaN(scanmin) || isNaN(scanmax)) {
    alert('Scan range is not numeric');
    return false;
  }
  if (isNaN(offmin) || isNaN(offmax)) {
    alert('Offset range is not numeric');
    return false;
  }
  if (document.trdsform.scan.value != 4) {
    if (document.trdsform.scan.value == 1) {
      if (scanmin < 0 || scanmin > 20 || scanmax < 0 || scanmax > 20) {
        alert('In this scanning mode scan angle must be in range [0:20] degr.');
        return false;
      }
      if (offmin < 0 || offmin > 20 || offmax < 0 || offmax > 20) {
        alert('In this scanning mode offset angle must be in range [0:20] degr.');
        return false;
      }
    }
    else if (document.trdsform.scan.value == 2) {
      if (scanmin < 0 || scanmin > 10 || scanmax < 0 || scanmax > 10) {
        alert('In this scanning mode scan angle must be in range [0:10] degr.');
        return false;
      }
      if (offmin < -10 || offmin > 10 || offmax < -10 || offmax > 10) {
        alert('In this scanning mode offset angle must be in range [-10:10] degr.');
        return false;
      }
    }
    else if (document.trdsform.scan.value == 3) {
      if (scanmin < 0 || scanmin > 20 || scanmax < 0 || scanmax > 20) {
        alert('In this scanning mode scan angle must be in range [0:20] degr.');
        return false;
      }
      if (offmin < 0 || offmin > 20 || offmax < 0 || offmax > 20) {
        alert('In this scanning mode offset angle must be in range [0:20] degr.');
        return false;
      }
    }
  } else {
    if (scanmin < -1 || scanmin > 1 || scanmax < -1 || scanmax > 1) {
      alert('In this scanning mode scan limits must be in range [-1:+1] A^-1.');
      return false;
    }
    if (offmin < 0 || offmin > 1 || offmax < 0 || offmax > 1) {
      alert('In this scanning mode offset limits must be in range [0:1] A^-1');
      return false;
    }
  }

  if (! valid_chars(document.trdsform.profile.value)) {
    alert('Illegal characters in the surface profile input.\n'+
          'Illegal characters are non-English letters, any\n' +
          'characters encoded in UTF-8, and "$\'<>|\\ \n' +
          'Possibly try to disable browser translating.');
    return false;
  }

  if (! valid_code_chars(document.trdsform.code.value)) {
    alert('Illegal characters in the crystal/material code input:\n\n'+
          document.trdsform.code.value + '\n\n' +
          'Illegal characters are non-English letters, any\n' +
          'characters encoded in UTF-8, and the special symbols\n' +
          'other than "()+-.". Possibly try to disable translating\n' +
          'into your language in the browser.');
    return false;
  }

// alert('OK '); return false;
   return true;
}
