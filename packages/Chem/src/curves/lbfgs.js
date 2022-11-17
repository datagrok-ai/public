function copyInto(target, source) {
	for (var i = 0; i < source.length; i++) {
		target[i] = source[i];
	}
}

function limitedMemoryBFGS(optimizable, parameters) {
	var lbfgsStart = +new Date();

	var converged = false;
	var maxIterations = 1000;
	var tolerance = 0.0001;
	var gradientTolerance = 0.001;
	var epsilon = 0.00001;
	
	var maxIterations = 100;
	var memorySize = 4;
	
	var numParameters = parameters.length;
	
	var gradient = numeric.rep([numParameters], 0.0);
	optimizable.getGradient(parameters, gradient);
	var oldGradient = numeric.clone(gradient);
	var oldParameters = numeric.clone(parameters);

	var direction = numeric.clone(gradient);
	
	// Project direction to the l2 ball
	numeric.diveq(direction, numeric.norm2(direction));
	
	var parameterChangeBuffer = []; // "s"
	var gradientChangeBuffer = []; // "y"
	var scaleBuffer = []; // "rho"
	
	// Initial step, do a line search in the direction of the gradient
	var scale = backtrackingLineSearch(optimizable, direction, gradient, parameters);
	
	// "parameters" has now been updated, so get a new value and gradient
	var value = optimizable.getValue(parameters);
	gradient = optimizable.getGradient(parameters, gradient);
	
	var oldValue = value;
	
	if (scale == 0.0) {
		console.log("Line search can't step in initial direction.");
	}
	
	for (var iteration = 0; iteration < maxIterations; iteration++) {
		var start = +new Date();
		var end;
		
		console.log("Beginning L-BFGS iteration, v=" + value + " ||g||=" + numeric.norm2(gradient));
		
		// Update the buffers with diffs
		if (parameterChangeBuffer.length < memorySize) {
			// If the buffer isn't full yet, add new arrays
			parameterChangeBuffer.unshift(numeric.sub(parameters, oldParameters));
			gradientChangeBuffer.unshift(numeric.sub(gradient, oldGradient));
		}
		else {
			// Otherwise, reuse the memory from the last array
			var parameterChange = parameterChangeBuffer.pop();
			var gradientChange = gradientChangeBuffer.pop();
			for (var i = 0; i < numParameters; i++) {
				parameterChange[i] = parameters[i] - oldParameters[i];
				gradientChange[i] = gradient[i] - oldGradient[i];
			}
			parameterChangeBuffer.unshift(parameterChange);
			gradientChangeBuffer.unshift(gradientChange);
		}
				
		// Save the old values. Gradient will be overwritten, then parameters.
		copyInto(oldParameters, parameters);
		copyInto(oldGradient, gradient);
	
		var sy = 0.0;
		var yy = 0.0;
		for (var i = 0; i < numParameters; i++) {
			sy += parameterChangeBuffer[0][i] * gradientChangeBuffer[0][i];
			yy += gradientChangeBuffer[0][i] * gradientChangeBuffer[0][i];
		}
		var scalingFactor = sy / yy;
		scaleBuffer.unshift(1.0 / sy);
		
		if (scalingFactor > 0.0) { console.log("Scaling factor greater than zero: " + scalingFactor); }
		
		// Renaming the "gradient" array to "direction" -- but it's the same memory.
		copyInto(direction, gradient);
		
		// Forward pass, from newest to oldest
		var alpha = [];
		for (var step = 0; step < parameterChangeBuffer.length; step++) {
			var currentAlpha = 0.0;
			for (var i = 0; i < numParameters; i++) {
				currentAlpha += parameterChangeBuffer[step][i] * direction[i];
			}
			currentAlpha *= scaleBuffer[step];
			
			//var currentAlpha = scaleBuffer[step] * numeric.dot(parameterChangeBuffer[step], direction)
			alpha.push(currentAlpha);
			for (var i = 0; i < numParameters; i++) {
				direction[i] += gradientChangeBuffer[step][i] * -currentAlpha;
			}
		}
		
		for (var i = 0; i < numParameters; i++) {
			direction[i] *= scalingFactor;
		}
		
		// Backward pass, from oldest to newest
		for (var step = parameterChangeBuffer.length - 1; step >= 0; step--) {
			var beta = 0.0;
			for (var i = 0; i < numParameters; i++) {
				beta += gradientChangeBuffer[step][i] * direction[i];
			}
			beta *= scaleBuffer[step];
			
			var currentAlpha = alpha[step];
			for (var i = 0; i < numParameters; i++) {
				direction[i] += parameterChangeBuffer[step][i] * (currentAlpha - beta);
			}
		}
		
		// Negate the direction, to maximize rather than minimize
		for (var i = 0; i < numParameters; i++) {
			direction[i] = -direction[i];
		}
		
		scale = backtrackingLineSearch(optimizable, direction, gradient, parameters);
		if (scale == 0.0) {
			console.log("Cannot step in current direction");
		}
		
		value = optimizable.getValue(parameters);
		gradient = optimizable.getGradient(parameters, gradient);
		
		// Test for convergence
		if (2.0 * (value - oldValue) <= tolerance * (Math.abs(value) + Math.abs(oldValue) + epsilon)) {
			console.log("Value difference below threshold: " + value + " - " + oldValue);
			end = +new Date();
			console.log("Finished iterations " + (end - lbfgsStart));
			return true;
		}
		
		var gradientNorm = numeric.norm2(gradient);
		if (gradientNorm < gradientTolerance) {
			console.log("Gradient norm below threshold: " + gradientNorm);
			end = +new Date();
			console.log("Finished iterations " + (end - lbfgsStart));
			return true;
		}
		else if (gradientNorm == 0.0) {
			console.log("Gradient norm is zero");
			end = +new Date();
			console.log("Finished iterations " + (end - lbfgsStart));
			return true;
		}
		
		oldValue = value;
		
	}

	end = +new Date();
	console.log("Finished iterations " + (end - lbfgsStart));
	return true;
}

function backtrackingLineSearch(optimizable, direction, gradient, parameters) {
	var numParameters = parameters.length;
	
	var MAXIMUM_STEP = 100.0;
	var RELATIVE_TOLERANCE = 0.0001;
	var DECREASE_FRACTION = 0.0001;
		
	var oldScale = 0.0;
	var scale = 1.0;
	var newScale = 0.0;
	
	var originalValue = optimizable.getValue(parameters);
	var oldValue = originalValue;
	
	// Make sure the initial step size isn't too big
	var twoNorm = numeric.norm2(direction);
	if (twoNorm > MAXIMUM_STEP) {
		console.log("Initial step " + twoNorm + " is too big, reducing")
		numeric.muleq(direction, MAXIMUM_STEP / twoNorm);
	}
	
	// Get the initial slope of the function of the scale.
	var slope = 0.0;
	for (var i = 0; i < numParameters; i++) {
		slope += gradient[i] * direction[i];
	}

	// Find the minimum acceptable scale value.
	var maxValue = 0.0;
	for (var i = 0; i < numParameters; i++) {
		var v = Math.abs( direction[i] / Math.max(Math.abs(parameters[i]), 1.0) );
		if (v > maxValue) { maxValue = v; }
	}
	var minimumScale = RELATIVE_TOLERANCE / maxValue;
	
	for (var iteration = 0; iteration < 25; iteration++) {
		for (var i = 0; i < numParameters; i++) {
			parameters[i] += (scale - oldScale) * direction[i];
		}
		
		if (scale < minimumScale) {
			console.log("Step too small, exiting.");
			return 0.0;
		}
		
		var value = optimizable.getValue(parameters);
		
		if (value >= originalValue + DECREASE_FRACTION * scale * slope) {
			//console.log("Exiting line search at value " + value);
			return scale;
		}
		else if (! isFinite(value)) {
			newScale = 0.2 * scale;
		}
		else {
			if (scale == 1.0) {
				// This is only true if this is the first iteration (?)
				newScale = -slope / (2.0 * (value - originalValue - slope));
			}
			else {
				var x1 = value - originalValue - scale * slope;
				var x2 = oldValue - originalValue - oldScale * slope;
				var oneOverScaleSquared = 1.0 / (scale * scale);
				var oneOverOldScaleSquared = 1.0 / (oldScale * oldScale);
				var oneOverScaleDiff = 1.0 / (scale - oldScale);
				
				var a = oneOverScaleDiff * (x1 * oneOverScaleSquared - x2 * oneOverOldScaleSquared);
				var b = oneOverScaleDiff * (-x1 * oldScale * oneOverScaleSquared + x2 * scale * oneOverOldScaleSquared);
				
				if (a == 0.0) {
					newScale = -slope / (2.0 * b);
				}
				else {
					var disc = b * b - 3.0 * a * slope;
					if (disc < 0.0) { newScale = 0.5 * scale; }
					else if (b <= 0.0) { newScale = (-b + Math.sqrt(disc)) / (3.0 * a); }
					else { newScale = -slope / (b + Math.sqrt(disc)); }
				}
				
				if (newScale > 0.5 * scale) { newScale = 0.5 * scale; }
			}
			
		}
		
		oldValue = value;
		oldScale = scale;
		scale = Math.max(newScale, 0.1 * scale);
	}
}

var quadratic = {
	getValue: function(parameters) {
		var x = parameters[0];
		var y = parameters[1];

		return -3*x*x - 4*y*y + 2*x - 4*y + 18;
	},

	getGradient: function (parameters, gradient) {
		gradient[0] = -6 * parameters[0] + 2;
		gradient[1] = -8 * parameters[1] - 4;
		return gradient;
    }
};

function doubleExp (n) {
	var x = Array(n);
	for (var i = 0; i < n; i++) {
		x[i] = Math.log(Math.random()) * ( Math.random() > 0.5 ? 1.0 : -1.0 );
	}
	return x;
}

var ridgeRegression = {
	covariates: [],
	responses: [],
	originalParameters: doubleExp(100),
	
	sample: function (n, noise) {
		for (var i = 0; i < n; i++) {
			var x = doubleExp(100);
			this.responses.push(numeric.dot(x, this.originalParameters) + noise());
			this.covariates.push(x);
		}
	},
	
	getValue: function(parameters) {
		var logLikelihood = 0.0;
		for (var i = 0; i < this.covariates.length; i++) {
			var residual = this.responses[i] - numeric.dot(this.covariates[i], parameters);
			logLikelihood += -0.5 * residual * residual;
		}
		return logLikelihood;
	},
	
	getGradient: function(parameters, gradient) {
		for (var i = 0; i < this.covariates.length; i++) {
			var residual = this.responses[i] - numeric.dot(this.covariates[i], parameters);
			numeric.addeq(gradient, numeric.mul(this.covariates[i], residual));
		}
		return gradient;
	}
};