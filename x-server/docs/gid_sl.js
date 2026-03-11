// GID_form fields enabling & disabling + fields validation
// Version-2.0 2013/05/19  fixed setOnloads as it did not work properly in IE after hitting "back" button

function LargeFloat(webpage, windowname, pagewidth, pageheight) {
  launch = window.open (webpage, windowname, 'width='+pagewidth+',height='+pageheight
                       +',toolbar=0,location=0,directories=0,status=0,menubar=0,scrollbars=1,resizable=1');
}

function setOnloads() {
  switchXrayFields();
  switchFcentre();
  switchScanAxis();
  switchStandingWaves();
}

function switchXrayFields() {
  if (document.gidform.xway.value != 3) {
     document.gidform.wave.disabled=false;
     document.gidform.line.disabled=true;
  } else {
     document.gidform.wave.disabled=true;
     document.gidform.line.disabled=false;
  }
}

function switchFcentre() {
  if (document.gidform.igie.value == 3 ||	//[3]. Surface orientation & condition of coplanar grazing incidence
      document.gidform.igie.value == 4 ||	//[4]. Surface orientation & condition of coplanar grazing exit
      document.gidform.igie.value == 5) {	//[5]. Surface orientation & condition of symmetric Bragg case
     // fcentre is a geometry parameter: [1,7]=incidence angle, [2,8]=exit angle, [6]=Bragg planes angle, [9]=g0/gh
     document.gidform.fcentre.disabled=true;
     document.gidform.unic.disabled=true;
  } else {
     document.gidform.fcentre.disabled=false;
     if (document.gidform.igie.value == 9) {document.gidform.unic.disabled=true;}
     else                                  {document.gidform.unic.disabled=false;}
  }
  if (document.gidform.igie.value <= 5) {
     document.gidform.n1.disabled=false;
     document.gidform.n2.disabled=false;
     document.gidform.n3.disabled=false;
     document.gidform.m1.disabled=false;
     document.gidform.m2.disabled=false;
     document.gidform.m3.disabled=false;
     document.gidform.miscut.disabled=false;
     document.gidform.unim.disabled=false;
  } else {
     document.gidform.n1.disabled=true;
     document.gidform.n2.disabled=true;
     document.gidform.n3.disabled=true;
     document.gidform.m1.disabled=true;
     document.gidform.m2.disabled=true;
     document.gidform.m3.disabled=true;
     document.gidform.miscut.disabled=true;
     document.gidform.unim.disabled=true;
  }
}

function switchScanAxis() {
  if (document.gidform.axis.value == 5) {
     document.gidform.a1.disabled=false;
     document.gidform.a2.disabled=false;
     document.gidform.a3.disabled=false;
  } else {
     document.gidform.a1.disabled=true;
     document.gidform.a2.disabled=true;
     document.gidform.a3.disabled=true;
  }
  if (document.gidform.axis.value >= 7) {
     document.gidform.unis.value=5;
     document.gidform.unis.disabled=true;
  } else {
     if (document.gidform.unis.value == 5) {document.gidform.unis.value=3;}
     document.gidform.unis.disabled=false;
  }
}

function switchStandingWaves() {
  if (document.gidform.swflag !== undefined) {
     if (document.gidform.swflag.checked) {
        document.gidform.swref.disabled=false;
        document.gidform.swmin.disabled=false;
        document.gidform.swmax.disabled=false;
        document.gidform.swpts.disabled=false;
        document.gidform.swphase.disabled=false;
     } else {
        document.gidform.swref.disabled=true;
        document.gidform.swmin.disabled=true;
        document.gidform.swmax.disabled=true;
        document.gidform.swpts.disabled=true;
        document.gidform.swphase.disabled=true;
     }
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

function gid_validate() {
  var wave_to_energy = 12.3981;
  var wave_min = 0.1;
  var wave_max = 10.;
  var energy_min = wave_to_energy / wave_max;
  var energy_max = wave_to_energy / wave_min;
  if (document.gidform.xway.value == 1) {				//Wavelength (A)
    if (document.gidform.wave.disabled) {
      // document.gidform.wave.disabled=false;
      alert('X-ray wavelength input is disabled by the browser Javascript. \n'+
            'This is unexpected behaviour (possibly untested browser?). \n'+
            'Please report the problem along with the version \n'+
            'of your web browser to X-ray Server administrator.');
      return false;
    }
    var wave = parseFloat(document.gidform.wave.value);
    if (wave < wave_min || wave > wave_max) {
      alert('X-ray wavelength='+wave+' Angstrom must be in range '+wave_min+'-'+wave_max+' Angstrom');
      return false;
    }
  }
  else if (document.gidform.xway.value == 2) {				//Energy (KeV)
    if (document.gidform.wave.disabled) {
      // document.gidform.wave.disabled=false;
       alert('X-ray energy input is disabled by the browser Javascript. \n'+
             'This is unexpected behaviour (possibly untested browser?). \n'+
             'Please report the problem along with the version \n'+
             'of your web browser to X-ray Server administrator.');
      return false;
    }
    var energy = parseFloat(document.gidform.wave.value);
    if (energy < energy_min || energy > energy_max) {
      alert('X-ray energy='+energy+' KeV must be in range: '
           +energy_min.toFixed(4)+'-'+energy_max.toFixed(2)+' KeV \n'
           +'(wavelength: '+wave_min+'-'+wave_max+' Angstrom)');
      return false;
    }
  }
  else if (document.gidform.xway.value == 3) {				//X-ray line
    if (document.gidform.line.disabled) {
      // document.gidform.line.disabled=false;
       alert('X-ray line selection is disabled by the browser Javascript. \n'+
             'This is unexpected behaviour (possibly untested browser?). \n'+
             'Please report the problem along with the version \n'+
             'of your web browser to X-ray Server administrator.');
      return false;
    }
    var selectedLine = document.gidform.line.selectedIndex;
    if (document.gidform.line.options[selectedLine].text == '') {
      alert('X-ray line selected, but not chosen from the list');
      return false;
    }
  }
  else if (document.gidform.xway.value == 4) {				//Bragg angle (degr)
    if (document.gidform.wave.disabled) {
      // document.gidform.wave.disabled=false;
       alert('X-ray Bragg angle input is disabled by the browser Javascript. \n'+
             'This is unexpected behaviour (possibly untested browser?). \n'+
             'Please report the problem along with the version \n'+
             'of your web browser to X-ray Server administrator.');
      return false;
    }
    if (document.gidform.wave.value <= 0. || document.gidform.wave.value > 90.) {
      alert('Bragg angle must be in range: 0-90 degr.');
      return false;
    }
  }
  if (document.gidform.i1.value.length == 0 ||
      document.gidform.i2.value.length == 0 ||
      document.gidform.i3.value.length == 0) {
    alert('Please fill non-zero Bragg reflection indices');
    return false;
  }
  if (document.gidform.i1.value == 0 &&
      document.gidform.i2.value == 0 &&
      document.gidform.i3.value == 0) {
    alert('Please specify non-zero Bragg reflection indices');
    return false;
  }
  if (document.gidform.w0.value <= 0 || document.gidform.w0.value > 99) {
    alert('The Debye-Waller modifier of x0 should be in range: [0-99], w0>0');
    return false;
  }
  if (document.gidform.wh.value <= 0 || document.gidform.wh.value > 99) {
    alert('The Debye-Waller modifier of xh should be in range: [0-99], wh>0');
    return false;
  }
  if (document.gidform.daa.value < -.5 || document.gidform.daa.value > .5) {
    alert('The da/a parameter of substrate should be in range: |da/a|<0.5');
    return false;
  }
  var fcentre = 0;
  if (! document.gidform.fcentre.disabled) {fcentre = parseFloat(document.gidform.fcentre.value);}
  if (document.gidform.igie.value != 3 &&		// coplanar grazing incidence(fcentre not used)
      document.gidform.igie.value != 4 &&		// coplanar grazing exit(fcentre not used)
      document.gidform.igie.value != 5 &&		// symmetric Bragg case (fcentre not used)
      document.gidform.igie.value != 9) {		// coplanar via beta (no units for fcentre)
    if (fcentre != 0 && document.gidform.unic.value == -1) {
      alert('Please specify units for geometry parameter');
      return false;
    }
    if (document.gidform.igie.value == 1 ||		// non-coplanar case via the incidence angle of k0
        document.gidform.igie.value == 2 ||		// non-coplanar case via the exit angle of kh
        document.gidform.igie.value == 6 ||		// coplanar case specified via the angle of Bragg planes to the surface
        document.gidform.igie.value == 7 ||		// coplanar case specified via the incidence angle of k0
        document.gidform.igie.value == 8) {		// coplanar case specified via the exit angle of kh
      var cfactor = 1.;
      if      (document.gidform.unic.value == 0) {cfactor = 1.;}		//degr.
      else if (document.gidform.unic.value == 1) {cfactor = 60.;}		//min.
      else if (document.gidform.unic.value == 2) {cfactor = 17.4532925;}	//mrad
      else if (document.gidform.unic.value == 3) {cfactor = 3600;}		//sec
      else if (document.gidform.unic.value == 4) {cfactor = 17453.2925;}	//urad
      fcentre = fcentre / cfactor;
    }
  }
  if (! document.gidform.miscut.disabled) {
    if (document.gidform.miscut.value.length == 0) {
      alert('Please fill the surface miscut angle');
      return false;
    }
    var mfactor = 1.;
    if      (document.gidform.unim.value == 0) {mfactor = 1.;}		//degr.
    else if (document.gidform.unim.value == 1) {mfactor = 60.;}		//min.
    else if (document.gidform.unim.value == 2) {mfactor = 17.4532925;}	//mrad
    else if (document.gidform.unim.value == 3) {mfactor = 3600;}	//sec
    else if (document.gidform.unim.value == 4) {mfactor = 17453.2925;}	//urad
    var miscut = parseFloat(document.gidform.miscut.value) / mfactor;
    if (isNaN(miscut)) {
      alert('Surface miscut angle is not numeric');
      return false;
    }
    if (miscut < -90 || miscut > 90) {
      alert('Surface miscut angle must be in range [-90:90] degr.');
      return false;
    }
    if (document.gidform.n1.value.length == 0 ||
        document.gidform.n2.value.length == 0 ||
        document.gidform.n3.value.length == 0) {
      alert('Please fill non-zero base surface normal indices');
      return false;
    }
    if (document.gidform.n1.value == 0 &&
        document.gidform.n2.value == 0 &&
        document.gidform.n3.value == 0) {
      alert('Please specify non-zero base surface normal indices');
      return false;
    }
    if (document.gidform.m1.value.length == 0 ||
        document.gidform.m2.value.length == 0 ||
        document.gidform.m3.value.length == 0) {
      alert('Please fill non-zero surface miscut direction indices');
      return false;
    }
    if (document.gidform.m1.value == 0 &&
        document.gidform.m2.value == 0 &&
        document.gidform.m3.value == 0) {
      alert('Please specify non-zero surface miscut direction indices');
      return false;
    }
    var nlen2 = document.gidform.n1.value * document.gidform.n1.value
              + document.gidform.n2.value * document.gidform.n2.value
              + document.gidform.n3.value * document.gidform.n3.value;
    var mlen2 = document.gidform.m1.value * document.gidform.m1.value
              + document.gidform.m2.value * document.gidform.m2.value
              + document.gidform.m3.value * document.gidform.m3.value;
    var nmcos = document.gidform.n1.value * document.gidform.m1.value
              + document.gidform.n2.value * document.gidform.m2.value
              + document.gidform.n3.value * document.gidform.m3.value;
    nmcos = nmcos * nmcos / (nlen2 * mlen2);
    if (Math.abs(nmcos) > 0.99) {
      alert('Surface miscut direction is parallel to base normal. Please correct.');
      return false;
    }
  }

  if (document.gidform.igie.value == 1) {		// non-coplanar via angle of k0
    if (! document.gidform.fcentre.disabled && document.gidform.fcentre.value.length == 0) {
      alert('Please fill the angle of k0 to the surface');
      return false;
    }
    else if (fcentre < -90 || fcentre > 90) {
      alert('Non-coplanar geometry specification: Incidence angle of k0 must be in range [-90:90] degr.');
      return false;
    }
  }
  else if (document.gidform.igie.value == 2) {		// non-coplanar via angle of kh
    if (! document.gidform.fcentre.disabled && document.gidform.fcentre.value.length == 0) {
      alert('Please fill the angle of kh to the surface');
      return false;
    }
    else if (fcentre < -90 || fcentre > 90) {
      alert('Non-coplanar geometry specification: Exit angle of kh must be in range [-90:90] degr.');
      return false;
    }
  }
  else if (document.gidform.igie.value == 6) {		// coplanar via Bragg planes angle
    if (! document.gidform.fcentre.disabled && document.gidform.fcentre.value.length == 0) {
      alert('Please fill the angle of Bragg planes to the surface');
      return false;
    }
    else if (fcentre < -90 || fcentre > 90) {
      alert('Coplanar geometry specification: the angle of Bragg planes to the surface must be in range [-90:90] degr.');
      return false;
    }
  }
  else if (document.gidform.igie.value == 7) {		// coplanar via incidence angle
    if (! document.gidform.fcentre.disabled && document.gidform.fcentre.value.length == 0) {
      alert('Please fill the incidence angle');
      return false;
    }
    else if (fcentre < -90 || fcentre > 90) {
      alert('Coplanar geometry specification: Incidence angle must be in range [-90:90] degr.');
      return false;
    }
  }
  else if (document.gidform.igie.value == 8) {		// coplanar via exit angle
    if (! document.gidform.fcentre.disabled && document.gidform.fcentre.value.length == 0) {
      alert('Please fill the exit angle');
      return false;
    }
    else if (fcentre < -90 || fcentre > 90) {
      alert('Coplanar geometry specification: Exit angle must be in range [-90:90] degr.');
      return false;
    }
  }
  else if (document.gidform.igie.value == 9) {
    if (! document.gidform.fcentre.disabled && document.gidform.fcentre.value.length == 0) {
      alert('Please fill the asymmetry factor, beta=g0/|gh|');
      return false;
    }
    if (fcentre < 0 || fcentre > 1.E5) {
      alert('Asymmetry factor should be in range 0-1E5');
      return false;
    }
  }
  if (document.gidform.axis.value == 5 &&
      document.gidform.a1.value == 0 &&
      document.gidform.a2.value == 0 &&
      document.gidform.a3.value == 0) {
    alert('Other scan axis is selected, but axis indices are not specified');
    return false;
  }
  var factor = 1.;
  if (document.gidform.axis.value <= 6) {				//angle scanning
    if      (document.gidform.unis.value == 0) {factor = 1.;}		//degr.
    else if (document.gidform.unis.value == 1) {factor = 60.;}		//min.
    else if (document.gidform.unis.value == 2) {factor = 17.4532925;}	//mrad
    else if (document.gidform.unis.value == 3) {factor = 3600;}		//sec
    else if (document.gidform.unis.value == 4) {factor = 17453.2925;}	//urad
  }
  var scanmin = parseFloat(document.gidform.scanmin.value) / factor;
  var scanmax = parseFloat(document.gidform.scanmax.value) / factor;
  if (isNaN(scanmin)) {
    alert('Minimum scan range is not numeric');
    return false;
  }
  if (isNaN(scanmax)) {
    alert('Maximum scan range is not numeric');
    return false;
  }
  if (document.gidform.axis.value <= 6) {
    var Scan_Max = 45;							//in degr.
    if (scanmin < -Scan_Max || scanmin > Scan_Max ||
        scanmax < -Scan_Max || scanmax > Scan_Max) {
      alert('Scan angle must be in range [-'+Scan_Max+':'+Scan_Max+'] degr.');
      return false;
    }
  } else {
    var Scan_Max = 1000;						//in eV
    if (scanmin < -Scan_Max || scanmin > Scan_Max ||
        scanmax < -Scan_Max || scanmax > Scan_Max) {
      alert('Scan energy must be in range [-'+Scan_Max+':'+Scan_Max+'] eV');
      return false;
    }
  }
  var nscanMax = 9999;
  if (document.gidform.nscan.value < 1 || document.gidform.nscan.value > nscanMax) {
    alert('Number of scan points must be in range [1-'+nscanMax+']');
    return false;
  }
  var alphamax_min=10;
  if (document.gidform.alphamax.value < alphamax_min && document.gidform.alphamax.value != 0) {
    alert('Alpha_max should be > '+alphamax_min);
    return false;
  }

  if (! valid_chars(document.gidform.profile.value)) {
    alert('Illegal characters in the surface profile input.\n'+
          'Illegal characters are non-English letters, any\n' +
          'characters encoded in UTF-8, and "$\'<>|\\ \n' +
          'Possibly try to disable browser translating.');
    return false;
  }

  if (! valid_code_chars(document.gidform.code.value)) {
    alert('Illegal characters in the crystal/material code input:\n\n'+
          document.gidform.code.value + '\n\n' +
          'Illegal characters are non-English letters, any\n' +
          'characters encoded in UTF-8, and the special symbols\n' +
          'other than "()+-.". Possibly try to disable translating\n' +
          'into your language in the browser.');
    return false;
  }

// if (document.gidform.codes_cryst.value != null) {document.gidform.codes_cryst.value = null;}
// if (document.gidform.codes_amo.value   != null) {document.gidform.codes_amo.value   = null;}
// if (document.gidform.codes_atom.value  != null) {document.gidform.codes_atom.value  = null;}

// alert('OK  ' + document.gidform.nscan.value); return false;
   return true;
}
