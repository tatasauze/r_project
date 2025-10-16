% SIDINIT  - Inicializa las variables globales del paquete de subspace identification.
% Debe ejecutarse una vez al principio de la sesi�n
% 23/12/96
% (C@@)
global SIDOPTION SIDERROR SIDWARN SIDDISP
if isempty(SIDOPTION)
   SIDOPTION = zeros(1,10);
end
sidopt;

% Inicializaci�n de los mensajes de error
SIDERROR = ['1. Ejecute SIDINIT antes de usar la librer�a                         ';
            '2. N�mero de argumentos incorrecto                                   ';
            '3. SIDOPT: Opci�n no reconocida: %s                                  ';
            '4. SIDOPT: Valor no reconocido: %s                                   ';
            '5. n debe ser menor o igual a i y a j                                '];
            
% Inicializaci�n de los mensajes de aviso
SIDWARN = [];

% Inicializaci�n de mensajes de informaci�n
SIDDISP = [];

