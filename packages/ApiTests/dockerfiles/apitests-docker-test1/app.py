import os

from flask import Flask, request

app = Flask('apitests')


@app.route('/square', methods=['GET'])
def square():
    num = int(request.args.get("number", 5))
    return {"result": num * num}, 200


@app.route('/health', methods=['GET'])
def health():
    num = 2 * 2
    if num == 4:
        return {"result": "Healthy"}, 200
    else:
        return {"result": "Unhealthy"}, 500


if __name__ == '__main__':
    if os.name == 'nt':
        temp_dir = os.environ.get('Temp')
    else:
        temp_dir = '/tmp'
    pid = open(os.path.join(temp_dir, 'apitests_pid.txt'), 'w')
    pid.write(str(os.getpid()))
    pid.close()

    app.run(host='0.0.0.0', port=5353, threaded=True)
