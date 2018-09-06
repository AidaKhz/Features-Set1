function yp = decimateFD(y, m) 
%   DECIMATE (DOWNSAMPLE) A SIGNAL IN FREQUENCY DOMAIN
%
%   Frequency domain decimation function to reduce the original sampling
%   rate of a signal to a lower rate.
%
%   Syntax:    
%       yp = decimateFD(y, m)
%
%   Input: 
%         y = input signal to be decimated (y must be a row vector)
%
%         m = decimation factor (> 1.0) (e.g., 2)
%
%        yp = decimated signal
%
%   Example: Input is a combined sinusoidal signal. It will
%   be downsampled by 2 (this example is modified from MatLAB)
%
%       m = 2;
%       t = 0:.00025:1;
%       x = sin(2*pi*30*t) + sin(2*pi*60*t);
%       subplot 211
%       stem(0:120,x(1:121),'filled','markersize',3)
%       xlabel 'Sample number',ylabel 'Original'; grid on
%       subplot 212
%       yp = decimateFD(x,m);
%       stem(0:60,yp(1:61),'filled','markersize',3,'color','r')
%       xlabel 'Sample number',ylabel 'Decimated'; grid on; 
%       print('-dpng','-painters','-r600','decimateFD.png')
%
%   See also resampleFDZP (frequency-domain zero padding upsampling)   
%
%   Note: No anti-aliasing filter is implemented. 
%
%   THIS SOFTWARE IS PROVIDED "AS IS" AND ANY EXPRESS OR IMPLIED
%   WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
%   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN
%   NO EVENT SHALL THE COPYRIGHT OWNER BE LIABLE FOR ANY DIRECT, INDIRECT,
%   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
%   OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
%   TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%   USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%   DAMAGE.
%
%   Written by Dr. Erol Kalkan, P.E. (ekalkan@usgs.gov)
%   $Revision: 1.0 $  $Date: 2017/05/31 14:03:00 $

% compute spectrum
z = fftshift(fft(y))/numel(y);

% half of total number of reduced samples
n = numel(z)*(m-1)/2/m; 

% drop n number of samples from beginning and end of spectrum 
z(1:n) = [];
z(end-n+1:end) = [];

% Compute inverse FFT by forcing conjugate symmetric
yp = ifft(ifftshift(z)*numel(z),'symmetric');